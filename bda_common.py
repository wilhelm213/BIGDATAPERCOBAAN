"""
bda_common.py — fondasi bersama untuk seluruh notebook pipeline BDA Influenza A.

Dipakai oleh 01_ingestion sampai 06_mart_dashboard. Setiap notebook memanggil
modul ini, jadi konfigurasi, path, sesi Spark, audit, dan monitoring hanya
didefinisikan satu kali.

Pemakaian di notebook:
    import sys; sys.path.insert(0, r"D:\\BDA\\nb")
    from bda_common import *
"""

from __future__ import annotations

import gc
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import time
import warnings
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

warnings.filterwarnings("ignore")

# Konsol Windows memakai cp1252. Data NCBI mengandung nama negara, galur, dan
# penyetor beraksen, sehingga print bisa gagal dengan UnicodeEncodeError di
# terminal maupun nbconvert. Paksa UTF-8 sekali di sini untuk seluruh notebook.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# ══════════════════════════════════════════════════════════════════════
# MODE — lokal (Windows) atau klaster (di dalam container Docker)
# ══════════════════════════════════════════════════════════════════════
# Notebook yang sama dipakai di kedua mode. Yang berbeda hanya tempat data
# disimpan dan ke mana Spark terhubung, dan itu ditentukan oleh variabel
# lingkungan yang disetel docker-compose.
# .strip() pada setiap nilai bukan kehati-hatian berlebihan. Di cmd.exe,
# `set NAMA=nilai && perintah` menyimpan "nilai " LENGKAP DENGAN spasi
# sebelum &&. Spasi itu tidak terlihat di mana pun, tetapi membuat Spark
# menolak master dengan pesan "Master must either be yarn or start with
# spark, mesos, k8s, or local" -- yang lalu muncul di notebook sebagai
# JAVA_GATEWAY_EXITED, galat yang sama sekali tidak menyebut penyebabnya.
MODE = os.environ.get("BDA_MODE", "lokal").strip().lower()
KLASTER = MODE == "klaster"
HDFS_URI = os.environ.get("BDA_HDFS", "hdfs://namenode:8020").strip()

# ══════════════════════════════════════════════════════════════════════
# PATH — zona data lake
# ══════════════════════════════════════════════════════════════════════
BASE = Path("/workspace") if KLASTER else Path(r"D:\BDA")

LAKE = BASE / "lake"        # zona bronze: hasil ingestion apa adanya
STAGE = BASE / "stage"      # zona silver: bersih, terstruktur, Parquet
FEAT = BASE / "features"    # feature store k-mer
MODELS = BASE / "models"    # model terlatih + metrik
GRAPH = BASE / "graph"      # simpul, tepi, hasil algoritma graf
MART = BASE / "mart"        # zona gold: tabel siap dashboard
OUT = BASE / "output"       # gambar, tabel, laporan

TMP = BASE / "spark_tmp"
CKPT = BASE / "checkpoint"
AUDIT = BASE / "audit"
LOGDIR = BASE / "spark_events"

for _d in (LAKE, STAGE, FEAT, MODELS, GRAPH, MART, OUT, TMP, CKPT, AUDIT, LOGDIR):
    _d.mkdir(parents=True, exist_ok=True)

# subfolder zona bronze, satu per jalur ingestion
LAKE_FASTA = LAKE / "fasta"          # jalur A — sekuens tak terstruktur
LAKE_META = LAKE / "meta_csv"        # jalur B — metadata terstruktur
for _d in (LAKE_FASTA, LAKE_META):
    _d.mkdir(parents=True, exist_ok=True)


def jalur(zona: str) -> str:
    """Alamat sebuah zona untuk dibaca atau ditulis Spark.

    Mode klaster mengembalikan URI HDFS; mode lokal mengembalikan folder biasa.
    Notebook memakai fungsi ini alih-alih menuliskan path secara langsung,
    sehingga berpindah mode tidak menuntut penyuntingan satu sel pun.

        jalur("stage/sequences")
          klaster -> hdfs://namenode:8020/bda/stage/sequences
          lokal   -> D:\\BDA\\stage\\sequences
    """
    zona = zona.strip("/")
    if KLASTER:
        return f"{HDFS_URI}/bda/{zona}"
    return str(BASE / Path(zona))


# ══════════════════════════════════════════════════════════════════════
# KONFIGURASI
# ══════════════════════════════════════════════════════════════════════
# Disetel untuk mesin ini: i7-14700F 20 core / 28 thread, RAM 63,7 GB,
# NVMe Gen4. Sisakan 4 core dan ~24 GB untuk Windows, Jupyter, dan browser.
CFG = {
    "taxid": 11320,                 # Influenza A virus

    # ---------- Spark ----------
    # Mode lokal memakai satu JVM di Windows; mode klaster memakai tiga worker
    # Docker dengan 4 core dan 7 GB masing-masing.
    "spark_master": (os.environ.get("BDA_SPARK_MASTER",
                                    "spark://spark-master:7077").strip()
                     if KLASTER else "local[16]"),
    "driver_memory": "8g" if KLASTER else "32g",

    # 4g, bukan 6g. Spark meminta executor_memory DITAMBAH overhead sekitar
    # 10% (minimal 384 MB), dan totalnya harus muat di bawah
    # yarn.scheduler.maximum-allocation-mb. Dengan plafon 6144 MB, nilai 6g
    # meminta 6144 + 614 = 6758 MB dan ditolak dengan pesan
    # "Required executor memory ... is above the max threshold of this
    # cluster". 4096 + 410 = 4506 MB memberi ruang yang lapang.
    #
    # Plafonnya sendiri diturunkan dari 7168 ke 6144 untuk memberi tempat
    # bagi Neo4j di dalam batas RAM WSL 48 GB.
    "executor_memory": "4g",
    "executor_cores": 4,
    "shuffle_part": 256,
    "max_result": "4g",

    # ---------- ingestion ----------
    "efetch_batch": 500,            # sekuens per permintaan efetch
    "vv_chunk": 100_000,            # baris per potongan CSV vvsearch2
    "ncbi_delay": 0.34,             # tanpa API key: maksimum 3 permintaan/detik
    "ncbi_api_key": None,           # isi lewat bda_secret.json -> 10 permintaan/detik
    "http_timeout": 600,
    "http_retry": 6,

    # ---------- kualitas data ----------
    "min_len": 500,                 # 6,1% sekuens di bawah ini, terlalu pendek untuk k-mer
    "max_len": 2600,                # segmen influenza terpanjang 2.341 bp
    "max_ambigu_pct": 5.0,          # persentase maksimum huruf non-ACGT

    # ---------- fitur k-mer ----------
    "k_utama": 8,
    "k_ablasi": [4, 6, 8, 10],
    "vocab_max": 1 << 18,           # dipakai HashingTF saat 4**k terlalu besar

    # ---------- machine learning ----------
    "seed": 42,
    "n_subtipe_teratas": 20,        # sisanya digabung ke kelas "LAINNYA"
    "fraksi_uji": 0.2,
    "fraksi_skalabilitas": [0.01, 0.05, 0.25, 0.50, 1.00],

    # ---------- graph ----------
    "lsh_hash": 5,

    # Ambang 0,35 dipakai pada versi awal dan DUA KALI membuat disk penuh
    # sampai nol. Jarak Jaccard 0,35 berarti "simpan tiap pasang yang
    # berbagi 65% k-mer". Lintas seluruh arsip itu menyaring banyak, tetapi
    # graf dibangun PER SEGMEN -- dan di dalam satu segmen sekuens influenza
    # memang nyaris identik. Ribuan H3N2 setahun berbagi jauh di atas 65%,
    # sehingga ambang itu praktis tidak menyaring apa pun dan keluaran join
    # mendekati O(n^2).
    "lsh_ambang_jarak": 0.10,       # jarak Jaccard; makin kecil makin ketat

    # Tiga penahan tambahan. Yang paling ampuh adalah membatasi SIMPUL,
    # karena ongkos join tumbuh kuadratik terhadapnya; membatasi tepi
    # setelah tepinya terlanjur dibuat tidak menolong apa-apa.
    "lsh_simpul_maks": 15_000,      # sampel bila satu segmen lebih besar
    "lsh_tetangga_maks": 20,        # tepi maksimum per simpul
    "lsh_tepi_maks": 3_000_000,     # katup pengaman per segmen
}

# Segmen influenza A. Nomor segmen adalah label kelas untuk klasifikasi Tugas 1.
SEGMEN = {
    1: ("PB2", "polymerase basic 2"),
    2: ("PB1", "polymerase basic 1"),
    3: ("PA", "polymerase acidic"),
    4: ("HA", "hemagglutinin"),
    5: ("NP", "nucleoprotein"),
    6: ("NA", "neuraminidase"),
    7: ("M", "matrix"),
    8: ("NS", "nonstructural"),
}
SEG_EKSTERNAL = [4, 6]              # HA dan NA — penentu subtipe
SEG_INTERNAL = [1, 2, 3, 5, 7, 8]   # gen internal — dipakai model pembanding

# Pemetaan nilai kolom segmen yang tidak seragam di GenBank ke nomor 1-8.
PETA_SEGMEN = {
    "1": 1, "2": 2, "3": 3, "4": 4, "5": 5, "6": 6, "7": 7, "8": 8,
    "PB2": 1, "PB1": 2, "PA": 3, "HA": 4, "NP": 5, "NA": 6,
    "M": 7, "MA": 7, "M1": 7, "M2": 7, "MP": 7,
    "NS": 8, "NS1": 8, "NS2": 8,
    "SEGMENT 1": 1, "SEGMENT 2": 2, "SEGMENT 3": 3, "SEGMENT 4": 4,
    "SEGMENT 5": 5, "SEGMENT 6": 6, "SEGMENT 7": 7, "SEGMENT 8": 8,
    "RNA 1": 1, "RNA 2": 2, "RNA 3": 3, "RNA 4": 4,
    "RNA 5": 5, "RNA 6": 6, "RNA 7": 7, "RNA 8": 8,
}

RE_SUBTIPE = re.compile(r"\bH(\d{1,2})N(\d{1,2})\b", re.I)
RE_SEG_DEFLINE = re.compile(r"\bsegment\s+([1-8])\b", re.I)
RE_STRAIN = re.compile(r"\(?(A/[^()]{3,90}?)\s*(?:\((H\d+N\d+)\))?\)?[,\s]", re.I)


def normalisasi_subtipe(nilai) -> str | None:
    """Kembalikan bentuk baku H#N# atau None bila tidak dikenali."""
    if not nilai:
        return None
    m = RE_SUBTIPE.search(str(nilai))
    return f"H{int(m.group(1))}N{int(m.group(2))}" if m else None


def normalisasi_segmen(nilai) -> int | None:
    """Petakan nilai kolom segmen yang beragam ke nomor 1-8."""
    if nilai is None:
        return None
    s = str(nilai).strip().upper()
    if not s:
        return None
    return PETA_SEGMEN.get(s)


# ══════════════════════════════════════════════════════════════════════
# RAHASIA — kredensial MySQL dan API key NCBI
# ══════════════════════════════════════════════════════════════════════
SECRET_PATH = BASE / "nb" / "bda_secret.json"
SECRET_CONTOH = {
    "mysql": {
        "host": "127.0.0.1",
        "port": 3306,
        "user": "root",
        "password": "ISI_PASSWORD_MYSQL_ANDA",
        "database": "bda_influenza",
    },
    "ncbi_api_key": "",
}


def muat_secret() -> dict:
    """Baca kredensial. Bila belum ada, tulis berkas contoh dan beri instruksi."""
    if not SECRET_PATH.exists():
        SECRET_PATH.write_text(
            json.dumps(SECRET_CONTOH, indent=2), encoding="utf-8")
        raise FileNotFoundError(
            f"Berkas kredensial dibuat di {SECRET_PATH}.\n"
            "Isi password MySQL Anda, lalu jalankan ulang sel ini.\n"
            "API key NCBI opsional (gratis di https://www.ncbi.nlm.nih.gov/account/) "
            "dan menaikkan batas dari 3 menjadi 10 permintaan per detik.")
    s = json.loads(SECRET_PATH.read_text(encoding="utf-8"))
    if s.get("ncbi_api_key"):
        CFG["ncbi_api_key"] = s["ncbi_api_key"]
        CFG["ncbi_delay"] = 0.11

    # MySQL berjalan di Windows, bukan di dalam klaster. Dari sisi container,
    # 127.0.0.1 menunjuk ke container itu sendiri, jadi alamatnya harus
    # ditukar ke nama host khusus yang disediakan Docker Desktop.
    if KLASTER and s.get("mysql", {}).get("host") in ("127.0.0.1", "localhost"):
        s["mysql"] = dict(s["mysql"], host="host.docker.internal")
    return s


def mysql_url(secret: dict | None = None) -> str:
    s = (secret or muat_secret())["mysql"]
    return (f"mysql+pymysql://{s['user']}:{s['password']}"
            f"@{s['host']}:{s['port']}/{s['database']}?charset=utf8mb4")


def mysql_jdbc(secret: dict | None = None) -> tuple[str, dict]:
    """URL JDBC dan properti untuk Spark write."""
    s = (secret or muat_secret())["mysql"]
    url = (f"jdbc:mysql://{s['host']}:{s['port']}/{s['database']}"
           "?useUnicode=true&characterEncoding=UTF-8&rewriteBatchedStatements=true")
    return url, {"user": s["user"], "password": s["password"],
                 "driver": "com.mysql.cj.jdbc.Driver"}


# ══════════════════════════════════════════════════════════════════════
# MONITORING & AUDIT — lapisan yang melintasi seluruh pipeline
# ══════════════════════════════════════════════════════════════════════
RUN_ID = datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
_JEJAK: list[dict] = []


def _ram():
    try:
        import psutil
        v = psutil.virtual_memory()
        return v.used / 1e9, v.available / 1e9
    except Exception:
        return float("nan"), float("nan")


@contextmanager
def Tahap(nama: str, lapisan: str = "-"):
    """Ukur durasi dan pemakaian RAM satu tahap, lalu catat ke audit trail."""
    t0 = time.time()
    u0, _ = _ram()
    print("\n" + "-" * 68)
    print(f"[>] {lapisan} | {nama}")
    status = "OK"
    try:
        yield
    except Exception as ex:
        status = f"GAGAL {type(ex).__name__}"
        raise
    finally:
        dt = time.time() - t0
        u1, bebas = _ram()
        baris = {"run_id": RUN_ID, "lapisan": lapisan, "tahap": nama,
                 "status": status, "detik": round(dt, 2),
                 "ram_delta_gb": round(u1 - u0, 2),
                 "ram_bebas_gb": round(bebas, 1),
                 "waktu": datetime.now(timezone.utc).isoformat()}
        _JEJAK.append(baris)
        with open(AUDIT / f"{RUN_ID}.jsonl", "a", encoding="utf-8") as f:
            f.write(json.dumps(baris) + "\n")
        print(f"[<] {status} | {dt:,.1f} detik | RAM bebas {bebas:,.1f} GB")


def jejak_df():
    import pandas as pd
    return pd.DataFrame(_JEJAK)


# ══════════════════════════════════════════════════════════════════════
# SPARK
# ══════════════════════════════════════════════════════════════════════
_spark = None
_HADOOP_HOME = None


def path_pendek(p) -> str:
    """Bentuk 8.3 dari sebuah path Windows, mis. C:\\PROGRA~1.

    Ini bukan kosmetik. Peluncur Spark di Windows tidak mengutip path yang
    mengandung spasi, sehingga JVM mati sebelum melapor bila JAVA_HOME atau
    SPARK_HOME memuat spasi -- persis kasus "C:\\Users\\Willy Boen" dan
    "C:\\Program Files". Bentuk 8.3 tidak pernah berspasi, jadi masalahnya
    hilang tanpa perlu memindahkan apa pun.
    """
    p = str(p)
    if os.name != "nt":
        return p
    try:
        import ctypes
        buf = ctypes.create_unicode_buffer(1024)
        n = ctypes.windll.kernel32.GetShortPathNameW(p, buf, 1024)
        return buf.value if n else p
    except Exception:
        return p


def siapkan_java() -> str | None:
    """PySpark di Windows menolak jalan bila JAVA_HOME kosong atau berspasi.

    Di mesin ini JDK terpasang tetapi variabel lingkungannya tidak pernah
    disetel, sehingga gateway JVM mati sebelum sempat melapor. Fungsi ini
    mencarinya sendiri lalu menyetelnya dalam bentuk 8.3.
    """
    jh = os.environ.get("JAVA_HOME", "").strip()
    if jh and (Path(jh) / "bin" / "java.exe").exists():
        jh = path_pendek(jh)
        os.environ["JAVA_HOME"] = jh
        return jh

    kandidat = []
    for akar in (r"C:\Program Files\Java", r"C:\Program Files\Eclipse Adoptium",
                 r"C:\Program Files\Microsoft", r"C:\Program Files\Amazon Corretto",
                 r"C:\Program Files\Zulu"):
        p = Path(akar)
        if p.is_dir():
            kandidat += [d for d in p.iterdir()
                         if d.is_dir() and (d / "bin" / "java.exe").exists()]

    if not kandidat:
        print("  [!] JDK tidak ditemukan. Pasang JDK 11 atau 17, "
              "lalu setel JAVA_HOME.")
        return None

    def skor(d: Path):
        # Spark 3.5 resmi mendukung Java 8, 11, dan 17. Utamakan 17, lalu 11.
        nama = d.name.lower()
        for v, s in (("17", 3), ("11", 2), ("1.8", 1), ("-8", 1)):
            if v in nama:
                return s
        return 0

    pilih = path_pendek(sorted(kandidat, key=skor, reverse=True)[0])
    os.environ["JAVA_HOME"] = pilih
    os.environ["PATH"] = pilih + r"\bin" + os.pathsep + os.environ.get("PATH", "")
    print(f"  JAVA_HOME  : {pilih}")
    return pilih


def siapkan_spark_home() -> str:
    """Setel SPARK_HOME ke bentuk 8.3 dari paket pyspark yang terpasang.

    Tanpa ini, path instalasi yang memuat spasi membuat peluncuran gagal
    dengan pesan JAVA_GATEWAY_EXITED yang menyesatkan.
    """
    import pyspark
    sh = path_pendek(Path(pyspark.__file__).parent)
    os.environ["SPARK_HOME"] = sh
    print(f"  SPARK_HOME : {sh}")
    return sh


def siapkan_hadoop() -> str | None:
    """Sediakan HADOOP_HOME bila winutils.exe tersedia.

    Spark di Windows mencari winutils.exe untuk operasi berkas lokal
    tertentu. Tanpa itu Spark tetap berjalan tetapi membanjiri log dengan
    peringatan; bila folder yang disiapkan ada, dipakai.
    """
    hh = os.environ.get("HADOOP_HOME", "").strip()
    calon = [Path(hh)] if hh else []
    calon.append(BASE / "hadoop")
    for c in calon:
        if (c / "bin" / "winutils.exe").exists():
            os.environ["HADOOP_HOME"] = str(c)
            os.environ["hadoop.home.dir"] = str(c)
            os.environ["PATH"] = str(c / "bin") + os.pathsep + os.environ.get("PATH", "")
            return str(c)
    return None


def spark_session(nama: str, memori: str | None = None,
                  master: str | None = None, paket: list[str] | None = None,
                  konfig: dict | None = None):
    """Buat SparkSession. Satu notebook satu sesi; matikan sebelum notebook berikutnya.

    `konfig` menambahkan atau menimpa setelan Spark apa pun. Dipakai mode
    YARN, yang butuh setelan yang tidak relevan di mode Standalone --
    jumlah executor, memori ApplicationMaster, dan nama host driver agar
    container di NodeManager bisa menghubunginya balik.
    """
    global _spark
    if _spark is not None:
        return _spark

    global _HADOOP_HOME
    if not KLASTER:
        # Penyesuaian khusus Windows. Di dalam container semuanya sudah benar.
        siapkan_java()
        siapkan_spark_home()
        _HADOOP_HOME = siapkan_hadoop()

    from pyspark.sql import SparkSession

    # appName ikut ke baris perintah cmd.exe saat JVM diluncurkan, jadi
    # karakter seperti | & < > akan merusak peluncuran dengan pesan
    # JAVA_GATEWAY_EXITED yang tidak menjelaskan apa-apa. Bersihkan dulu.
    nama_aman = re.sub(r"[^A-Za-z0-9._-]+", "-", f"BDA-Influenza-{nama}").strip("-")

    b = (SparkSession.builder
         .appName(nama_aman)
         .master(master or CFG["spark_master"])
         .config("spark.driver.memory", memori or CFG["driver_memory"])
         .config("spark.driver.maxResultSize", CFG["max_result"])
         .config("spark.sql.shuffle.partitions", CFG["shuffle_part"])
         .config("spark.sql.adaptive.enabled", "true")
         .config("spark.sql.adaptive.coalescePartitions.enabled", "true")
         .config("spark.serializer", "org.apache.spark.serializer.KryoSerializer")
         .config("spark.sql.execution.arrow.pyspark.enabled", "true")
         .config("spark.ui.showConsoleProgress", "true"))

    if KLASTER:
        b = (b.config("spark.executor.memory", CFG["executor_memory"])
              .config("spark.executor.cores", str(CFG["executor_cores"]))
              # dfs.replication adalah setelan SISI KLIEN. Menyetelnya di
              # namenode saja tidak cukup: tanpa baris ini Spark menulis
              # dengan default 3 salinan, dan konsumsi ruang naik 50%.
              .config("spark.hadoop.dfs.replication", "2")
              .config("spark.eventLog.enabled", "true")
              .config("spark.eventLog.dir", f"{HDFS_URI}/bda/spark-events")
              # Driver berada di container jupyter; executor harus bisa
              # menghubunginya balik lewat nama host itu.
              .config("spark.driver.host", os.environ.get("HOSTNAME", "jupyter")))
    else:
        b = b.config("spark.local.dir", str(TMP))
        # Event log menulis lewat Hadoop FileSystem, dan di Windows itu
        # menuntut winutils.exe. Tanpa winutils, mengaktifkannya membuat
        # SparkContext gagal start -- jadi hanya dinyalakan bila winutils ada.
        if _HADOOP_HOME:
            b = (b.config("spark.eventLog.enabled", "true")
                  .config("spark.eventLog.dir", LOGDIR.as_uri()))
        else:
            print("  Event log NONAKTIF (winutils.exe tidak ada di mode lokal).")

    if paket:
        b = b.config("spark.jars.packages", ",".join(paket))

    # Diterapkan paling akhir supaya pemanggil bisa menimpa setelan mana pun
    # di atas, bukan sekadar menambah yang belum ada.
    for k, v in (konfig or {}).items():
        b = b.config(k, str(v))

    # Pemeriksaan antrean hanya berlaku untuk Spark Standalone, yang memberi
    # aplikasi pertama seluruh core. YARN punya penjadwal sendiri dan tidak
    # berperilaku begitu, jadi pemeriksaannya dilewati.
    if KLASTER and not str(master or CFG["spark_master"]).startswith("yarn"):
        periksa_klaster_kosong()

    _spark = b.getOrCreate()
    _spark.sparkContext.setLogLevel("WARN")
    print(f"Spark {_spark.version} | master {_spark.sparkContext.master} | "
          f"driver {memori or CFG['driver_memory']} | "
          f"paralelisme {_spark.sparkContext.defaultParallelism}")
    print(f"Spark UI: {_spark.sparkContext.uiWebUrl}")
    return _spark


def unggah_ke_hdfs(spark, sumber: Path, tujuan: str, pola: str = "*") -> int:
    """Salin berkas dari filesystem lokal container ke HDFS.

    Image jupyter tidak memuat perkakas baris perintah `hdfs`, tetapi Spark
    sudah membawa pustaka Hadoop lengkap. Jadi penyalinan dilakukan lewat
    API FileSystem di JVM, tanpa perlu memasang apa pun.

    Berkas yang sudah ada di tujuan dilewati, sehingga fungsi ini aman
    dijalankan berulang kali.
    """
    jvm = spark._jvm
    konf = spark._jsc.hadoopConfiguration()
    JPath = jvm.org.apache.hadoop.fs.Path
    fs = jvm.org.apache.hadoop.fs.FileSystem.get(jvm.java.net.URI(tujuan), konf)
    fs.mkdirs(JPath(tujuan))

    baru = 0
    berkas = sorted(p for p in sumber.glob(pola) if p.is_file())
    for f in berkas:
        dst = JPath(f"{tujuan}/{f.name}")
        if fs.exists(dst):
            continue
        # (delSrc=False, overwrite=True)
        fs.copyFromLocalFile(False, True, JPath(str(f)), dst)
        baru += 1
    print(f"  {sumber.name}: {baru} berkas baru diunggah, "
          f"{len(berkas) - baru} sudah ada -> {tujuan}")
    return baru


def periksa_klaster_kosong(diam: bool = False) -> bool:
    """Peringatkan bila klaster sudah dipakai aplikasi Spark lain.

    Spark standalone memberi aplikasi pertama SELURUH core yang tersedia
    secara bawaan. Kalau satu notebook masih memegang sesinya, notebook
    berikutnya tidak akan mendapat satu core pun -- ia hanya mencetak
    "Initial job has not accepted any resources" berulang kali tanpa pernah
    maju. Gejalanya mudah disalahartikan sebagai "lambat", padahal sebenarnya
    menunggu selamanya.

    Pemeriksaan ini menjadikannya jelas sejak awal.
    """
    import urllib.request
    try:
        with urllib.request.urlopen("http://spark-master:8080/json/", timeout=8) as r:
            d = json.load(r)
    except Exception:
        return True

    aktif = [a for a in d.get("activeapps", []) if a.get("cores", 0) > 0]
    bebas = d.get("cores", 0) - d.get("coresused", 0)
    if not aktif:
        return True

    if not diam:
        print("  [!] KLASTER SEDANG DIPAKAI aplikasi lain:")
        for a in aktif:
            print(f"        {a['name']}  core={a['cores']}  "
                  f"jalan {a['duration']/1000/60:.0f} menit")
        print(f"      Core bebas: {bebas} dari {d.get('cores', 0)}")
        if bebas == 0:
            print("      Sesi ini akan MENGANTRE tanpa batas waktu.")
            print("      Hentikan sesi notebook lain dulu (jalankan stop_spark()")
            print("      di notebook itu, atau Kernel -> Shutdown Kernel).")
    return bebas > 0


# ══════════════════════════════════════════════════════════════════════
# YARN — manajemen sumber daya klaster
# ══════════════════════════════════════════════════════════════════════
YARN_RM = os.environ.get("BDA_YARN_RM", "resourcemanager:8088")


def yarn_api(jalan: str, timeout: int = 8):
    """Ambil satu endpoint REST ResourceManager, kembalikan dict atau None.

    Dipakai untuk membuktikan bahwa job memang berjalan di atas YARN,
    bukan di mode `local` yang diam-diam menjadi bawaan MapReduce ketika
    mapreduce.framework.name tidak disetel.
    """
    import urllib.request
    try:
        with urllib.request.urlopen(
                f"http://{YARN_RM}/ws/v1/cluster/{jalan.strip('/')}",
                timeout=timeout) as r:
            return json.load(r)
    except Exception:
        return None


def info_yarn(diam: bool = False) -> dict | None:
    """Ringkas keadaan klaster YARN: node aktif, memori, dan aplikasi."""
    m = yarn_api("metrics")
    if m is None:
        if not diam:
            print(f"  [!] ResourceManager tidak terjangkau di {YARN_RM}.")
            print("      Klaster dinyalakan dengan profil yarn?")
            print("      docker compose --profile yarn up -d")
        return None

    d = m.get("clusterMetrics", {})
    if not diam:
        print(f"  ResourceManager : {YARN_RM}")
        print(f"  NodeManager     : {d.get('activeNodes', 0)} aktif, "
              f"{d.get('lostNodes', 0)} hilang")
        print(f"  Memori          : {d.get('totalMB', 0):,} MB total, "
              f"{d.get('availableMB', 0):,} MB bebas")
        print(f"  vCore           : {d.get('totalVirtualCores', 0)} total, "
              f"{d.get('availableVirtualCores', 0)} bebas")
        print(f"  Aplikasi        : {d.get('appsRunning', 0)} berjalan, "
              f"{d.get('appsCompleted', 0)} selesai, "
              f"{d.get('appsFailed', 0)} gagal")
    return d


def yarn_aplikasi(batas: int = 10):
    """Daftar aplikasi terakhir di YARN sebagai DataFrame.

    Inilah bukti paling langsung bahwa sebuah job benar-benar melewati
    YARN: kalau tidak muncul di sini, job itu tidak pernah menyentuh
    klaster sama sekali.
    """
    import pandas as pd
    d = yarn_api("apps")
    if not d or not d.get("apps"):
        return pd.DataFrame(columns=["id", "name", "applicationType",
                                     "state", "finalStatus"])
    baris = []
    for a in d["apps"].get("app", [])[:batas]:
        baris.append({
            "id": a.get("id"),
            "nama": a.get("name"),
            "jenis": a.get("applicationType"),
            "status": a.get("state"),
            "hasil": a.get("finalStatus"),
            "detik": round(a.get("elapsedTime", 0) / 1000, 1),
            "memori_MB_detik": a.get("memorySeconds"),
            "vcore_detik": a.get("vcoreSeconds"),
        })
    return pd.DataFrame(baris)


def stop_spark():
    """Matikan sesi dan JVM-nya. Panggil di sel terakhir setiap notebook."""
    global _spark
    if _spark is not None:
        try:
            _spark.stop()
        except Exception:
            pass
        _spark = None
    gc.collect()
    try:
        import psutil
        sisa = [p.pid for p in psutil.process_iter(["name"])
                if "java" in (p.info["name"] or "").lower()]
        print(f"Spark dihentikan. Proses Java tersisa: {len(sisa)}")
    except Exception:
        print("Spark dihentikan.")


def bersihkan_java():
    """Darurat: bunuh seluruh proses Java yang menggantung."""
    import subprocess
    import psutil
    pids = [p.pid for p in psutil.process_iter(["name"])
            if "java" in (p.info["name"] or "").lower()]
    if pids:
        subprocess.run(["taskkill", "/F", "/IM", "java.exe"],
                       capture_output=True, text=True)
        time.sleep(3)
    v = _ram()
    print(f"{len(pids)} proses Java dimatikan. RAM bebas {v[1]:,.1f} GB")


# ══════════════════════════════════════════════════════════════════════
# UTILITAS
# ══════════════════════════════════════════════════════════════════════
def sesi_http():
    import requests
    s = requests.Session()
    s.headers.update({"User-Agent": "BDA-Influenza-Pipeline/2.0 (akademik)"})
    return s


def ukuran(p: Path) -> str:
    """Ukuran berkas atau folder dalam satuan yang mudah dibaca."""
    if not p.exists():
        return "0 B"
    n = (p.stat().st_size if p.is_file()
         else sum(f.stat().st_size for f in p.rglob("*") if f.is_file()))
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:,.1f} {u}"
        n /= 1024
    return f"{n:,.1f} PB"


def ringkas_zona():
    """Tampilkan isi tiap zona data lake beserta ukurannya."""
    import pandas as pd
    baris = []
    for nama, p in [("lake/fasta", LAKE_FASTA), ("lake/meta_csv", LAKE_META),
                    ("stage", STAGE), ("features", FEAT), ("models", MODELS),
                    ("graph", GRAPH), ("mart", MART), ("output", OUT)]:
        n = len(list(p.rglob("*"))) if p.exists() else 0
        baris.append({"zona": nama, "isi": n, "ukuran": ukuran(p)})
    return pd.DataFrame(baris)


def simpan_ckpt(nama: str, data: dict):
    (CKPT / f"{nama}.json").write_text(json.dumps(data, indent=2), encoding="utf-8")


def muat_ckpt(nama: str, bawaan: dict | None = None) -> dict:
    p = CKPT / f"{nama}.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            pass
    return dict(bawaan or {})


def sha1(teks: str) -> str:
    return hashlib.sha1(teks.encode("utf-8", "ignore")).hexdigest()


def info_mesin():
    import psutil
    v = psutil.virtual_memory()
    d = shutil.disk_usage(str(BASE))
    print(f"CPU logis      : {os.cpu_count()}")
    print(f"RAM total      : {v.total / 1e9:,.1f} GB   bebas {v.available / 1e9:,.1f} GB")
    print(f"Disk D: bebas  : {d.free / 1e9:,.1f} GB dari {d.total / 1e9:,.1f} GB")
    print(f"Python         : {sys.version.split()[0]}")
    print(f"Mode           : {MODE.upper()}" + (f"  ({HDFS_URI})" if KLASTER else "  (Windows lokal)"))
    print(f"run_id         : {RUN_ID}")


__all__ = [
    "MODE", "KLASTER", "HDFS_URI", "jalur",
    "BASE", "LAKE", "LAKE_FASTA", "LAKE_META", "STAGE", "FEAT", "MODELS",
    "GRAPH", "MART", "OUT", "TMP", "CKPT", "AUDIT", "LOGDIR",
    "CFG", "SEGMEN", "SEG_EKSTERNAL", "SEG_INTERNAL", "PETA_SEGMEN",
    "RE_SUBTIPE", "RE_SEG_DEFLINE", "RE_STRAIN",
    "normalisasi_subtipe", "normalisasi_segmen",
    "muat_secret", "mysql_url", "mysql_jdbc", "SECRET_PATH",
    "Tahap", "jejak_df", "RUN_ID",
    "spark_session", "stop_spark", "bersihkan_java", "unggah_ke_hdfs",
    "periksa_klaster_kosong",
    "YARN_RM", "yarn_api", "info_yarn", "yarn_aplikasi",
    "sesi_http", "ukuran", "ringkas_zona", "simpan_ckpt", "muat_ckpt",
    "sha1", "info_mesin", "siapkan_java", "siapkan_hadoop",
    "siapkan_spark_home", "path_pendek",
    "Path", "json", "re", "time", "gzip", "gc", "os", "sys", "datetime", "timezone",
]
