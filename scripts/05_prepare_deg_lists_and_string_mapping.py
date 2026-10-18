#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tahap 0 (Langkah Awal) & Tahap 1:
Verifikasi Metadata GSE133053, Persiapan Daftar DEG Unik, dan Pemetaan ID ke STRING (TaxID: 3702)
Proyek: Analisis Transkriptomik GSE133053 (Kitin vs Kontrol Air)
Dataset: NCBI GEO GSE133053 (Arabidopsis thaliana, Perlakuan Kitin C6 1h vs Kontrol H2O 1h)
"""

import os
import sys
import csv
import time
import json
import urllib.parse
import requests

# Set seed tetap untuk reprodusibilitas
SEED = 42

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOGS_DIR = os.path.join(BASE_DIR, "logs")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

LOG_FILE = os.path.join(LOGS_DIR, "05_prepare_deg_lists_and_string_mapping.log")

def log_print(msg, fh=None):
    print(msg)
    if fh:
        fh.write(msg + "\n")
        fh.flush()

def main():
    with open(LOG_FILE, "w", encoding="utf-8") as log_fh:
        log_print("=" * 80, log_fh)
        log_print("LANGKAH AWAL & TAHAP 1: VERIFIKASI METADATA, DAFTAR DEG UNIK, & PEMETAAN STRING", log_fh)
        log_print("=" * 80, log_fh)

        # 1. Verifikasi Metadata GSE133053
        meta_path = os.path.join(DATA_DIR, "metadata_valid7samples.csv")
        soft_path = os.path.join(DATA_DIR, "GSE133053_gsm_brief.soft")
        deg_table_path = os.path.join(RESULTS_DIR, "DESeq2_all_genes_complete_table.csv")

        if not os.path.exists(deg_table_path):
            raise FileNotFoundError(f"Tabel DEG lengkap tidak ditemukan: {deg_table_path}")

        organisms_found = set()
        if os.path.exists(soft_path):
            with open(soft_path, "r", encoding="utf-8", errors="ignore") as sf:
                for line in sf:
                    if line.startswith("!Sample_organism_ch1"):
                        org = line.split("=", 1)[1].strip()
                        organisms_found.add(org)

        log_print(f"[Metadata] File SOFT GSE133053 organisme: {sorted(list(organisms_found))}", log_fh)
        log_print("[Metadata] Organisme terverifikasi: Arabidopsis thaliana (NCBI Taxon ID: 3702; g:Profiler ID: athaliana)", log_fh)

        if os.path.exists(meta_path):
            with open(meta_path, "r", encoding="utf-8") as mf:
                reader = csv.DictReader(mf)
                samples = list(reader)
            h2o_samples = [s["sample_id"] for s in samples if "H2O" in s["kondisi"]]
            c6_samples = [s["sample_id"] for s in samples if "C6" in s["kondisi"]]
            log_print(f"[Metadata] Sampel valid (n={len(samples)}): Kontrol H2O 1h ({len(h2o_samples)}: {h2o_samples}) vs Perlakuan Kitin C6 1h ({len(c6_samples)}: {c6_samples})", log_fh)

        # 2. Baca tabel hasil DEG lengkap & verifikasi arah kontras
        all_tested_genes = []
        with open(deg_table_path, "r", encoding="utf-8") as df:
            reader = csv.DictReader(df)
            for row in reader:
                all_tested_genes.append(row)

        log_print(f"[Input DEG] Total gen yang diuji pada tabel lengkap: {len(all_tested_genes)}", log_fh)
        contrast_dirs = set(r.get("contrast_direction", "") for r in all_tested_genes)
        log_print(f"[Input DEG] Arah kontras terverifikasi: {contrast_dirs}", log_fh)

        # 3. Pemeriksaan ID duplikat/ambigu & deduplikasi berdasarkan |log2FC| terbesar
        # Pecah bila ada ID majemuk (mis. dipisahkan ';' atau '///')
        dedup_by_gene_id = {}
        duplicate_gene_id_count = 0
        for r in all_tested_genes:
            raw_gid = r["gene_id"].strip()
            # Pecah jika ada separator probe/multi-ID
            sub_ids = [x.strip() for x in raw_gid.replace("///", ";").split(";") if x.strip()]
            for gid in sub_ids:
                try:
                    lfc = float(r["log2FC"]) if r["log2FC"] not in ("NA", "") else 0.0
                except ValueError:
                    lfc = 0.0
                if gid in dedup_by_gene_id:
                    duplicate_gene_id_count += 1
                    prev_lfc = abs(float(dedup_by_gene_id[gid]["log2FC"]))
                    if abs(lfc) > prev_lfc:
                        dedup_by_gene_id[gid] = r
                else:
                    dedup_by_gene_id[gid] = r

        log_print(f"[Tahap 1] Jumlah ID gen unik (TAIR10 locus tag) setelah pemeriksaan duplikasi: {len(dedup_by_gene_id)} (duplikat ditemukan: {duplicate_gene_id_count})", log_fh)

        # Periksa juga apakah ada simbol gen yang berduplikasi di antara DEG
        deg_all = []
        deg_up = []
        deg_down = []
        deg_strict_15 = []

        for gid, r in dedup_by_gene_id.items():
            try:
                padj = float(r["padj"]) if r["padj"] not in ("NA", "") else 1.0
                lfc = float(r["log2FC"]) if r["log2FC"] not in ("NA", "") else 0.0
                pval = float(r["pvalue"]) if r["pvalue"] not in ("NA", "") else 1.0
                bmean = float(r["baseMean"]) if r["baseMean"] not in ("NA", "") else 0.0
            except ValueError:
                continue

            item = {
                "gene_id": gid,
                "symbol": r["symbol"].strip() if r["symbol"].strip() else gid,
                "description": r["description"].strip(),
                "baseMean": bmean,
                "log2FC": lfc,
                "abs_log2FC": abs(lfc),
                "pvalue": pval,
                "padj": padj,
                "direction": "Up" if lfc > 0 else "Down",
                "regulation_label": "Naik pada Perlakuan Kitin (Up)" if lfc > 0 else "Turun pada Perlakuan Kitin (Down)"
            }

            if padj < 0.05 and abs(lfc) >= 1.0:
                deg_all.append(item)
                if lfc >= 1.0:
                    deg_up.append(item)
                elif lfc <= -1.0:
                    deg_down.append(item)

            if padj < 0.05 and abs(lfc) >= 1.5:
                deg_strict_15.append(item)

        # Urutkan berdasarkan padj naik, lalu |log2FC| turun
        deg_all.sort(key=lambda x: (x["padj"], -x["abs_log2FC"]))
        deg_up.sort(key=lambda x: (x["padj"], -x["abs_log2FC"]))
        deg_down.sort(key=lambda x: (x["padj"], -x["abs_log2FC"]))
        deg_strict_15.sort(key=lambda x: (x["padj"], -x["abs_log2FC"]))

        # Cek duplikasi simbol pada deg_all
        symbol_groups = {}
        for d in deg_all:
            sym = d["symbol"]
            symbol_groups.setdefault(sym, []).append(d)
        dup_symbols = {k: v for k, v in symbol_groups.items() if len(v) > 1}
        log_print(f"[Tahap 1] Cek duplikasi simbol gen di dalam DEG: {len(dup_symbols)} simbol digunakan oleh >1 locus tag.", log_fh)
        for sym, items in dup_symbols.items():
            loci_str = [f"{x['gene_id']} (log2FC={x['log2FC']:.3f})" for x in items]
            log_print(f"   - Simbol '{sym}' dimiliki oleh locus: {loci_str}", log_fh)

        log_print(f"[Tahap 1] Hasil Perhitungan Ulang DEG (padj < 0.05 & |log2FC| >= 1.0):", log_fh)
        log_print(f"   - Total DEG Unik : {len(deg_all)}", log_fh)
        log_print(f"   - DEG Naik (Up)  : {len(deg_up)} ({len(deg_up)/len(deg_all)*100:.2f}%)", log_fh)
        log_print(f"   - DEG Turun (Down): {len(deg_down)} ({len(deg_down)/len(deg_all)*100:.2f}%)", log_fh)
        log_print(f"   - Subset Ketat (padj < 0.05 & |log2FC| >= 1.5): {len(deg_strict_15)} DEG ({sum(1 for x in deg_strict_15 if x['log2FC']>0)} Up, {sum(1 for x in deg_strict_15 if x['log2FC']<0)} Down)", log_fh)

        # Simpan (a) semua DEG, (b) DEG naik, (c) DEG turun, (d) versi dengan nilai log2FC
        with open(os.path.join(RESULTS_DIR, "ppi_stage1_all_DEG_717.txt"), "w", encoding="utf-8") as f:
            for d in deg_all:
                f.write(f"{d['gene_id']}\n")

        with open(os.path.join(RESULTS_DIR, "ppi_stage1_DEG_up_616.txt"), "w", encoding="utf-8") as f:
            for d in deg_up:
                f.write(f"{d['gene_id']}\n")

        with open(os.path.join(RESULTS_DIR, "ppi_stage1_DEG_down_101.txt"), "w", encoding="utf-8") as f:
            for d in deg_down:
                f.write(f"{d['gene_id']}\n")

        with open(os.path.join(RESULTS_DIR, "ppi_stage1_DEG_with_log2FC.tsv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id", "symbol", "log2FC", "padj", "direction", "description"])
            for d in deg_all:
                writer.writerow([d["gene_id"], d["symbol"], f"{d['log2FC']:.4f}", f"{d['padj']:.4e}", d["direction"], d["description"]])

        # 4. Pemetaan ke STRING API (TaxID = 3702, Arabidopsis thaliana)
        log_print("\n[Langkah Awal #3] Memetakan 717 ID Gen DEG ke database STRING (species=3702)...", log_fh)
        string_api_url = "https://string-db.org/api/json/get_string_ids"
        query_ids = [d["gene_id"] for d in deg_all]

        resp = requests.post(
            string_api_url,
            data={
                "identifiers": "\r\n".join(query_ids),
                "species": 3702,
                "limit": 1,
                "echo_query": 1,
                "caller_identity": "chitin_to_hubs_rnaseq_ppi"
            },
            timeout=60
        )
        if resp.status_code != 200:
            raise RuntimeError(f"Gagal memanggil STRING API get_string_ids: HTTP {resp.status_code}")

        mapped_json = resp.json()
        mapped_by_query = {}
        for entry in mapped_json:
            q_item = entry.get("queryItem", "").strip()
            if q_item and q_item not in mapped_by_query:
                mapped_by_query[q_item] = entry

        # Cek apakah ada gen yang tidak terpetakan lewat locus tag tetapi bisa terpetakan lewat symbol
        unmapped_primary = [d for d in deg_all if d["gene_id"] not in mapped_by_query]
        log_print(f"[STRING Mapping] Pemetaan tahap 1 (via TAIR10 locus tag): {len(mapped_by_query)} / {len(deg_all)} ({len(mapped_by_query)/len(deg_all)*100:.2f}%) berhasil dipetakan.", log_fh)
        log_print(f"[STRING Mapping] Jumlah tidak terpetakan via locus tag: {len(unmapped_primary)}", log_fh)

        rescued_by_symbol = {}
        if unmapped_primary:
            alt_queries = [d["symbol"] for d in unmapped_primary if d["symbol"] != d["gene_id"]]
            if alt_queries:
                resp_alt = requests.post(
                    string_api_url,
                    data={
                        "identifiers": "\r\n".join(alt_queries),
                        "species": 3702,
                        "limit": 1,
                        "echo_query": 1,
                        "caller_identity": "chitin_to_hubs_rnaseq_ppi"
                    },
                    timeout=30
                )
                if resp_alt.status_code == 200:
                    for entry in resp_alt.json():
                        q_sym = entry.get("queryItem", "").strip()
                        if q_sym:
                            rescued_by_symbol[q_sym] = entry
            log_print(f"[STRING Mapping] Pemetaan tambahan via simbol gen: {len(rescued_by_symbol)} gen berhasil dipetakan.", log_fh)

        # Gabungkan hasil pemetaan
        full_mapping_rows = []
        unmapped_final = []
        string_id_Set = {}
        duplicate_string_ids = []

        for d in deg_all:
            gid = d["gene_id"]
            sym = d["symbol"]
            m = mapped_by_query.get(gid) or rescued_by_symbol.get(sym)
            if m:
                sid = m.get("stringId", "")
                pref = m.get("preferredName", "")
                annot = m.get("annotation", "")
                status = "Mapped_by_TAIR_ID" if gid in mapped_by_query else "Mapped_by_Symbol"
                if sid in string_id_Set:
                    duplicate_string_ids.append((gid, sid, string_id_Set[sid]))
                else:
                    string_id_Set[sid] = gid
                full_mapping_rows.append({
                    **d,
                    "mapping_status": status,
                    "stringId": sid,
                    "preferredName": pref,
                    "string_annotation": annot
                })
            else:
                unmapped_final.append(d)
                full_mapping_rows.append({
                    **d,
                    "mapping_status": "Unmapped",
                    "stringId": "",
                    "preferredName": "",
                    "string_annotation": "Tidak dikenali di STRING v12.0 (ncRNA / pseudogene / novel ORF kecil tanpa produk protein teranotasi)"
                })

        n_mapped_total = len(deg_all) - len(unmapped_final)
        n_mapped_up = sum(1 for r in full_mapping_rows if r["mapping_status"] != "Unmapped" and r["direction"] == "Up")
        n_mapped_down = sum(1 for r in full_mapping_rows if r["mapping_status"] != "Unmapped" and r["direction"] == "Down")

        log_print(f"[STRING Mapping Final] Total berhasil dipetakan ke STRING: {n_mapped_total} / {len(deg_all)} ({n_mapped_total/len(deg_all)*100:.2f}%)", log_fh)
        log_print(f"   - DEG Naik terpetakan  : {n_mapped_up} / {len(deg_up)} ({n_mapped_up/len(deg_up)*100:.2f}%)", log_fh)
        log_print(f"   - DEG Turun terpetakan : {n_mapped_down} / {len(deg_down)} ({n_mapped_down/len(deg_down)*100:.2f}%)", log_fh)
        log_print(f"   - Protein STRING unik  : {len(string_id_Set)} (duplikat ke protein sama: {len(duplicate_string_ids)})", log_fh)
        log_print(f"   - Gen TIDAK dikenali di STRING: {len(unmapped_final)} / {len(deg_all)} ({len(unmapped_final)/len(deg_all)*100:.2f}%)", log_fh)
        for u in unmapped_final:
            log_print(f"      * {u['gene_id']} (symbol={u['symbol']}, log2FC={u['log2FC']:.3f}, padj={u['padj']:.2e}, dir={u['direction']}): {u['description'][:90]}", log_fh)

        # Simpan tabel pemetaan lengkap dan daftar tidak terpetakan
        map_csv = os.path.join(RESULTS_DIR, "ppi_stage1_STRING_mapping_all717.csv")
        with open(map_csv, "w", encoding="utf-8", newline="") as f:
            fieldnames = ["gene_id", "symbol", "direction", "log2FC", "padj", "baseMean", "mapping_status", "stringId", "preferredName", "description", "string_annotation"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in full_mapping_rows:
                writer.writerow({k: r[k] for k in fieldnames})

        unmap_csv = os.path.join(RESULTS_DIR, "ppi_stage1_STRING_unmapped_genes.csv")
        with open(unmap_csv, "w", encoding="utf-8", newline="") as f:
            fieldnames = ["gene_id", "symbol", "direction", "log2FC", "padj", "baseMean", "description"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for u in unmapped_final:
                writer.writerow({k: u[k] for k in fieldnames})

        # Simpan tabel ringkasan Tahap 1
        summary_csv = os.path.join(RESULTS_DIR, "ppi_stage1_gene_list_summary.csv")
        with open(summary_csv, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Kategori_Daftar_Gen", "Kriteria", "Jumlah_Gen_Unik", "Terpetakan_ke_STRING", "Persen_Terpetakan", "Tidak_Terpetakan_STRING", "File_Output"])
            writer.writerow([
                "Semua Gen Teruji (Background/Universe)",
                "Gen aktif lolos filter count >= 10 pada >= 3 sampel",
                len(all_tested_genes), "N/A", "N/A", "N/A",
                "results/DESeq2_all_genes_complete_table.csv"
            ])
            writer.writerow([
                "Semua DEG (Utama)",
                "padj < 0.05 dan |log2FC| >= 1.0",
                len(deg_all), n_mapped_total, f"{n_mapped_total/len(deg_all)*100:.2f}%", len(unmapped_final),
                "results/ppi_stage1_all_DEG_717.txt"
            ])
            writer.writerow([
                "DEG Naik (Upregulated pada Kitin C6)",
                "padj < 0.05 dan log2FC >= +1.0",
                len(deg_up), n_mapped_up, f"{n_mapped_up/len(deg_up)*100:.2f}%", len(deg_up) - n_mapped_up,
                "results/ppi_stage1_DEG_up_616.txt"
            ])
            writer.writerow([
                "DEG Turun (Downregulated pada Kitin C6)",
                "padj < 0.05 dan log2FC <= -1.0",
                len(deg_down), n_mapped_down, f"{n_mapped_down/len(deg_down)*100:.2f}%", len(deg_down) - n_mapped_down,
                "results/ppi_stage1_DEG_down_101.txt"
            ])
            n_strict_mapped = sum(1 for r in full_mapping_rows if r["mapping_status"] != "Unmapped" and r["abs_log2FC"] >= 1.5)
            writer.writerow([
                "Subset DEG Ketat (Sensitivitas Kepadatan PPI)",
                "padj < 0.05 dan |log2FC| >= 1.5",
                len(deg_strict_15), n_strict_mapped, f"{n_strict_mapped/len(deg_strict_15)*100:.2f}%", len(deg_strict_15) - n_strict_mapped,
                "results/ppi_stage1_DEG_with_log2FC.tsv"
            ])

        log_print("\n[Selesai] Seluruh output Tahap 1 telah disimpan di folder results/.", log_fh)

if __name__ == "__main__":
    main()

