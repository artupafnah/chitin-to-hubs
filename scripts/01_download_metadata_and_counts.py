#!/usr/bin/env python3
"""
Tahap 1 & 4: Pengunduhan Metadata SRA/GEO dan Matriks Raw Counts Resmi NCBI GEO
Dataset: GSE133053 / BioProject PRJNA549850 / SRA SRP201971
Publikasi: Wan et al. (2019) Nature Plants, PMID: 30252144
Organisme: Arabidopsis thaliana (ekotipe Col-0, kecambah / seedlings)
Kontras Utama: Perlakuan Elisitor Kitin (C6, 1 jam) vs Kontrol Air (H2O, 1 jam)
"""

import os
import re
import tarfile
import gzip
import urllib.request
import xml.etree.ElementTree as ET
import csv

# 1. Buat struktur direktori wajib sesuai aturan kerja
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIRS = ["data", "qc", "results", "figures", "scripts", "logs"]
for d in DIRS:
    os.makedirs(os.path.join(BASE_DIR, d), exist_ok=True)

print(f"[INFO] Direktori kerja: {BASE_DIR}")

# 2. Unduh metadata GSM dari NCBI GEO (SOFT format brief)
geo_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133053&targ=gsm&form=text&view=brief"
soft_path = os.path.join(BASE_DIR, "data", "GSE133053_gsm_brief.soft")
if not os.path.exists(soft_path):
    print("[INFO] Mengunduh metadata GSM dari NCBI GEO...")
    urllib.request.urlretrieve(geo_url, soft_path)

# Parse SOFT file
samples = []
current = {}
with open(soft_path, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.strip()
        if line.startswith("^SAMPLE = "):
            if current:
                samples.append(current)
            current = {"GSM": line.split("= ")[1].strip(), "chars": []}
        elif line.startswith("!Sample_title = "):
            current["title"] = line.split("= ", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = "):
            current["chars"].append(line.split("= ", 1)[1].strip())
        elif line.startswith("!Sample_instrument_model = "):
            current["instrument"] = line.split("= ", 1)[1].strip()
        elif line.startswith("!Sample_library_strategy = "):
            current["strategy"] = line.split("= ", 1)[1].strip()
        elif line.startswith("!Sample_relation = BioSample: "):
            current["BioSample"] = line.split("/")[-1].strip()
        elif line.startswith("!Sample_relation = SRA: "):
            current["SRX"] = line.split("term=")[-1].strip()
        elif line.startswith("!Sample_supplementary_file_1 = "):
            current["suppl_file"] = line.split("/")[-1].strip()
if current:
    samples.append(current)

for s in samples:
    for c in s.get("chars", []):
        if c.startswith("time point:"):
            s["time_point"] = c.split(":", 1)[1].strip()
        elif c.startswith("agent:"):
            s["agent"] = c.split(":", 1)[1].strip()
        elif c.startswith("ecotype:"):
            s["ecotype"] = c.split(":", 1)[1].strip()
    m = re.match(r"(\d+\s*h)\s+([A-Za-z0-9]+)_(\d+)", s.get("title", ""))
    if m:
        s["replicate"] = int(m.group(3))
    else:
        s["replicate"] = 1
    m_lane = re.search(r"RunId(\d+)_LaneId(\d+)", s.get("suppl_file", ""))
    if m_lane:
        s["run_id"] = m_lane.group(1)
        s["lane_id"] = f"Lane_{m_lane.group(2)}"
    else:
        s["run_id"] = "0195"
        s["lane_id"] = "Lane_1"

print(f"[INFO] Berhasil mem-parsing {len(samples)} sampel GSM dari GSE133053.")

# 3. Unduh SRA RunInfo dari NCBI Entrez untuk PRJNA549850 (SRP201971)
sra_runinfo_path = os.path.join(BASE_DIR, "data", "SraRunTable_PRJNA549850.csv")
if not os.path.exists(sra_runinfo_path):
    print("[INFO] Mengambil daftar ID SRA dari NCBI E-utilities (PRJNA549850)...")
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA549850[BioProject]&retmax=200"
    req = urllib.request.urlopen(esearch_url)
    xml_data = req.read()
    root = ET.fromstring(xml_data)
    id_list = [elem.text for elem in root.findall(".//Id")]
    print(f"[INFO] Ditemukan {len(id_list)} eksperimen/run SRA.")
    if id_list:
        ids_str = ",".join(id_list)
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={ids_str}&rettype=runinfo&retmode=text"
        urllib.request.urlretrieve(efetch_url, sra_runinfo_path)
        print(f"[INFO] SRA RunInfo disimpan ke {sra_runinfo_path}")

sra_by_srx = {}
with open(sra_runinfo_path, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("Run") and row.get("Run") != "Run":
            sra_by_srx[row.get("Experiment", "")] = row

print(f"[INFO] SRA RunInfo memuat {len(sra_by_srx)} baris SRR unik.")

# Gabungkan GSM dan SRA
merged_rows = []
for s in samples:
    srx = s.get("SRX", "")
    sra = sra_by_srx.get(srx, {})
    row = {
        "GSM": s.get("GSM", ""),
        "SRR": sra.get("Run", ""),
        "SRX": srx,
        "BioSample": s.get("BioSample", sra.get("BioSample", "")),
        "title": s.get("title", ""),
        "agent": s.get("agent", ""),
        "time_point": s.get("time_point", ""),
        "replicate": s.get("replicate", 1),
        "lane_id": s.get("lane_id", ""),
        "layout": sra.get("LibraryLayout", ""),
        "strategy": sra.get("LibraryStrategy", "RNA-Seq"),
        "selection": sra.get("LibrarySelection", "cDNA"),
        "instrument": sra.get("Model", s.get("instrument", "")),
        "spots_total_reads": int(sra.get("spots", 0)) if sra.get("spots", "").isdigit() else 0,
        "bases": int(sra.get("bases", 0)) if sra.get("bases", "").isdigit() else 0,
        "avgLength": sra.get("avgLength", ""),
        "size_MB": sra.get("size_MB", ""),
        "suppl_file": s.get("suppl_file", "")
    }
    merged_rows.append(row)

# Simpan semua 52 sampel
full_meta_path = os.path.join(BASE_DIR, "data", "metadata_GSE133053_all52samples.csv")
with open(full_meta_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(merged_rows[0].keys()))
    writer.writeheader()
    writer.writerows(merged_rows)
print(f"[INFO] Metadata lengkap 52 sampel disimpan di: {full_meta_path}")

# 4. Filter sampel untuk Kontras Utama:
# Perlakuan Kitin 1 jam (C6, 1 h) vs Kontrol Air 1 jam (H2O, 1 h)
focus_rows = [r for r in merged_rows if r["time_point"] == "1 h" and r["agent"] in ("C6", "H2O")]
focus_rows.sort(key=lambda x: (0 if x["agent"] == "H2O" else 1, int(x["replicate"])))

meta_out_rows = []
for r in focus_rows:
    kondisi = "Perlakuan_Kitin_C6" if r["agent"] == "C6" else "Kontrol_H2O"
    cond_short = "Chitin_1h" if r["agent"] == "C6" else "Control_1h"
    sample_id = f"{r['agent']}_1h_Rep{r['replicate']}"
    meta_out_rows.append({
        "sample_id": sample_id,
        "GSM": r["GSM"],
        "SRR": r["SRR"],
        "SRX": r["SRX"],
        "BioSample": r["BioSample"],
        "kondisi": kondisi,
        "condition_short": cond_short,
        "agent": r["agent"],
        "time_point": r["time_point"],
        "replicate": r["replicate"],
        "lane_id": r["lane_id"],
        "layout": r["layout"],
        "strand": "unstranded (Illumina TruSeq RNA v2)",
        "instrument": r["instrument"],
        "spots_total_reads": r["spots_total_reads"],
        "bases": r["bases"],
        "avgLength": r["avgLength"],
        "size_MB": r["size_MB"],
        "suppl_file": r["suppl_file"]
    })

meta_csv_path = os.path.join(BASE_DIR, "data", "metadata.csv")
with open(meta_csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(meta_out_rows[0].keys()))
    writer.writeheader()
    writer.writerows(meta_out_rows)

print("\n[INFO] Tabel metadata sampel utama (data/metadata.csv):")
for r in meta_out_rows:
    print(f"  {r['sample_id']:12s} | {r['GSM']} | {r['SRR']} | {r['kondisi']:18s} | Rep {r['replicate']} | {r['lane_id']} | {r['layout']} | Reads: {r['spots_total_reads']:,} | Len: {r['avgLength']} bp")

# 5. Verifikasi dan ekstrak GSE133053_RAW.tar lokal (TANPA unduh ulang)
raw_tar_path = os.path.join(BASE_DIR, "data", "GSE133053_RAW.tar")
raw_extract_dir = os.path.join(BASE_DIR, "data", "raw_counts_gsm")
os.makedirs(raw_extract_dir, exist_ok=True)

if not os.path.exists(raw_tar_path) or os.path.getsize(raw_tar_path) == 0:
    raise FileNotFoundError(f"File lokal {raw_tar_path} tidak ditemukan!")
print(f"\n[INFO] Menggunakan arsip lokal {raw_tar_path} ({os.path.getsize(raw_tar_path):,} bytes)")

with tarfile.open(raw_tar_path, "r") as tar:
    tar.extractall(path=raw_extract_dir)
print(f"[INFO] Berhasil mengekstrak file count per sampel ke {raw_extract_dir}")

# 6. Periksa struktur file count dan gabungkan menjadi matriks count (Tahap 4 Output)
sample_file = os.path.join(raw_extract_dir, meta_out_rows[0]["suppl_file"])
with gzip.open(sample_file, "rt") as gf:
    head_lines = [next(gf).strip() for _ in range(10)]
print("\n[INFO] Contoh 10 baris pertama file count:", meta_out_rows[0]["suppl_file"])
for hl in head_lines:
    print("  ", hl)
