#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tahap 3, 4, dan 5:
Jaringan Interaksi Protein-Protein (PPI) STRING, Analisis Centrality (NetworkX),
Deteksi Komunitas Louvain, Skor Prioritas Komposit, & Uji Stabilitas Peringkat
Proyek: Analisis Transkriptomik GSE133053 (Kitin vs Kontrol Air)
Dataset: NCBI GEO GSE133053 (Arabidopsis thaliana, TaxID: 3702, Kitin C6 1h vs Kontrol H2O 1h)
"""

import os
import sys
import csv
import math
import random
import requests
import numpy as np
import networkx as nx
from scipy.stats import spearmanr, kendalltau
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Set seed tetap untuk reprodusibilitas
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
LOGS_DIR = os.path.join(BASE_DIR, "logs")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

LOG_FILE = os.path.join(LOGS_DIR, "06_string_ppi_network_centrality_and_prioritization.log")

# Palet warna ramah buta warna (Okabe-Ito)
COLOR_UP = "#D55E00"      # Vermilion (Naik pada Kitin C6)
COLOR_DOWN = "#0072B2"    # Blue (Turun pada Kitin C6)
COLOR_TOP10_BORDER = "#F0E442" # Yellow/Gold highlight
COLOR_NEUTRAL = "#7F7F7F"

def log_print(msg, fh=None):
    print(msg)
    if fh:
        fh.write(msg + "\n")
        fh.flush()

def zscore_sample(arr):
    """Hitung Z-score menggunakan simpangan baku sampel (ddof=1, setara STDEV.S di Excel)."""
    a = np.asarray(arr, dtype=float)
    mean_val = np.mean(a)
    std_val = np.std(a, ddof=1)
    if std_val == 0 or np.isnan(std_val):
        return np.zeros_like(a)
    return (a - mean_val) / std_val

def compute_network_stats(G_active, total_deg_input, total_mapped_proteins, label):
    """Hitung statistik topologi jaringan lengkap."""
    all_mapped_nodes = total_mapped_proteins
    n_active = G_active.number_of_nodes()
    n_edges = G_active.number_of_edges()
    n_isolated = all_mapped_nodes - n_active

    avg_deg_all_mapped = (2.0 * n_edges / all_mapped_nodes) if all_mapped_nodes > 0 else 0.0
    avg_deg_active = (2.0 * n_edges / n_active) if n_active > 0 else 0.0
    density_active = nx.density(G_active) if n_active > 1 else 0.0

    if n_active > 0:
        comps = sorted(nx.connected_components(G_active), key=len, reverse=True)
        n_comps_active = len(comps)
        n_comps_total = n_comps_active + n_isolated
        gcc_nodes_set = comps[0]
        G_gcc = G_active.subgraph(gcc_nodes_set).copy()
        gcc_n = G_gcc.number_of_nodes()
        gcc_e = G_gcc.number_of_edges()
        gcc_frac = (gcc_n / n_active * 100.0) if n_active > 0 else 0.0
        gcc_clust = nx.average_clustering(G_gcc) if gcc_n > 1 else 0.0
        gcc_diam = nx.diameter(G_gcc) if gcc_n > 1 else 0
        gcc_aspl = nx.average_shortest_path_length(G_gcc) if gcc_n > 1 else 0.0
    else:
        n_comps_active = 0
        n_comps_total = n_isolated
        gcc_n = 0
        gcc_e = 0
        gcc_frac = 0.0
        gcc_clust = 0.0
        gcc_diam = 0
        gcc_aspl = 0.0
        G_gcc = nx.Graph()

    return {
        "Konfigurasi_Jaringan": label,
        "Total_DEG_Input": total_deg_input,
        "Protein_Terpetakan_STRING": all_mapped_nodes,
        "Node_Aktif_Berinteraksi": n_active,
        "Node_Terisolasi_Degree0": n_isolated,
        "Jumlah_Edge": n_edges,
        "Rata2_Degree_Semua_Protein_Terpetakan": round(avg_deg_all_mapped, 4),
        "Rata2_Degree_Node_Aktif": round(avg_deg_active, 4),
        "Kepadatan_Jaringan_Aktif": round(density_active, 5),
        "Jumlah_Komponen_Terhubung_Aktif": n_comps_active,
        "Jumlah_Komponen_Total_Termasuk_Isolat": n_comps_total,
        "Node_GCC_Komponen_Terbesar": gcc_n,
        "Edge_GCC_Komponen_Terbesar": gcc_e,
        "Persen_Node_Aktif_di_GCC": round(gcc_frac, 2),
        "Clustering_Coefficient_GCC": round(gcc_clust, 4),
        "Diameter_GCC": gcc_diam,
        "Rata2_Jarak_Terpendek_GCC": round(gcc_aspl, 4)
    }, G_gcc

def compute_centralities_and_composite(G, weighted=False):
    """
    Hitung 4 ukuran centrality (Degree, Betweenness, Closeness, Eigenvector) + Z-score (STDEV.S, ddof=1)
    + Composite Score (rata-rata 4 Z-score).
    Catatan metodologis:
    - Bila weighted=False: semua bobot dianggap 1 (topologi murni).
    - Bila weighted=True:
      * Degree (Strength) & Eigenvector menggunakan 'weight' = combined_score (kekuatan interaksi).
      * Betweenness & Closeness di NetworkX menginterpretasikan parameter weight/distance sebagai
        JARAK/BIAYA lintasan terpendek, sehingga digunakan 'distance' = 1.0 / combined_score.
    """
    nodes = sorted(list(G.nodes()))
    n = len(nodes)
    if n == 0:
        return {}

    if not weighted:
        deg_raw = {u: G.degree(u) for u in nodes}
        deg_cent = nx.degree_centrality(G)
        bet_cent = nx.betweenness_centrality(G, normalized=True, weight=None, seed=SEED)
        clo_cent = nx.closeness_centrality(G, distance=None, wf_improved=True)
        try:
            eig_cent = nx.eigenvector_centrality_numpy(G, weight=None)
        except Exception:
            eig_cent = nx.eigenvector_centrality(G, max_iter=2000, tol=1e-6, weight=None)
    else:
        deg_raw = {u: G.degree(u) for u in nodes}
        # Weighted degree centrality = sum(weight) / (N - 1)
        deg_cent = {u: sum(d.get("weight", 1.0) for _, _, d in G.edges(u, data=True)) / max(1, n - 1) for u in nodes}
        bet_cent = nx.betweenness_centrality(G, normalized=True, weight="distance", seed=SEED)
        clo_cent = nx.closeness_centrality(G, distance="distance", wf_improved=True)
        try:
            eig_cent = nx.eigenvector_centrality_numpy(G, weight="weight")
        except Exception:
            eig_cent = nx.eigenvector_centrality(G, max_iter=2000, tol=1e-6, weight="weight")

    arr_deg = np.array([deg_cent[u] for u in nodes])
    arr_bet = np.array([bet_cent[u] for u in nodes])
    arr_clo = np.array([clo_cent[u] for u in nodes])
    arr_eig = np.array([eig_cent[u] for u in nodes])

    z_deg = zscore_sample(arr_deg)
    z_bet = zscore_sample(arr_bet)
    z_clo = zscore_sample(arr_clo)
    z_eig = zscore_sample(arr_eig)

    comp_scores = (z_deg + z_bet + z_clo + z_eig) / 4.0

    results = {}
    for idx, u in enumerate(nodes):
        results[u] = {
            "degree_raw": deg_raw[u],
            "degree_centrality": float(arr_deg[idx]),
            "betweenness_centrality": float(arr_bet[idx]),
            "closeness_centrality": float(arr_clo[idx]),
            "eigenvector_centrality": float(arr_eig[idx]),
            "z_degree": float(z_deg[idx]),
            "z_betweenness": float(z_bet[idx]),
            "z_closeness": float(z_clo[idx]),
            "z_eigenvector": float(z_eig[idx]),
            "composite_score": float(comp_scores[idx])
        }

    # Urutkan menurun berdasarkan composite_score, tie-breaker degree_centrality
    ranked_nodes = sorted(nodes, key=lambda u: (results[u]["composite_score"], results[u]["degree_centrality"]), reverse=True)
    for r_idx, u in enumerate(ranked_nodes, start=1):
        results[u]["rank"] = r_idx

    return results

def save_fig_multi(fig, base_name):
    for ext in ["png", "pdf", "svg"]:
        out_path = os.path.join(FIGURES_DIR, f"{base_name}.{ext}")
        fig.savefig(out_path, dpi=300 if ext == "png" else None, bbox_inches="tight", facecolor="white")

def main():
    with open(LOG_FILE, "w", encoding="utf-8") as log_fh:
        log_print("=" * 85, log_fh)
        log_print("TAHAP 3, 4, & 5: JARINGAN PPI STRING, CENTRALITY, KOMUNITAS, & PRIORITISASI KANDIDAT", log_fh)
        log_print("=" * 85, log_fh)
        log_print(f"NetworkX version: {nx.__version__} | NumPy version: {np.__version__} | Seed: {SEED}", log_fh)

        # 1. Baca tabel pemetaan STRING dari Tahap 1
        map_csv = os.path.join(RESULTS_DIR, "ppi_stage1_STRING_mapping_all717.csv")
        if not os.path.exists(map_csv):
            raise FileNotFoundError(f"File pemetaan tidak ditemukan: {map_csv}")

        all_deg_rows = []
        mapped_rows = []
        string_to_gene = {}
        gene_meta = {}

        with open(map_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for r in reader:
                r["log2FC"] = float(r["log2FC"])
                r["padj"] = float(r["padj"])
                r["baseMean"] = float(r["baseMean"])
                r["abs_log2FC"] = abs(r["log2FC"])
                all_deg_rows.append(r)
                gene_meta[r["gene_id"]] = r
                if r["mapping_status"] != "Unmapped" and r["stringId"]:
                    mapped_rows.append(r)
                    string_to_gene[r["stringId"]] = r["gene_id"]

        log_print(f"[Input] Total DEG: {len(all_deg_rows)} | Terpetakan ke STRING: {len(mapped_rows)} protein unik.", log_fh)

        # 2. Panggil STRING API untuk mengambil interaksi jaringan pada skor >= 0.400 (mencakup >= 0.700)
        log_print("\n[Tahap 3] Mengambil jaringan interaksi PPI dari STRING API (species=3702, required_score=400)...", log_fh)
        net_url = "https://string-db.org/api/json/network"
        string_ids_list = [r["stringId"] for r in mapped_rows]

        resp = requests.post(
            net_url,
            data={
                "identifiers": "\r\n".join(string_ids_list),
                "species": 3702,
                "required_score": 400,
                "caller_identity": "chitin_to_hubs_rnaseq_ppi"
            },
            timeout=90
        )
        if resp.status_code != 200:
            raise RuntimeError(f"Gagal memanggil STRING API network: HTTP {resp.status_code}")

        raw_edges = resp.json()
        log_print(f"[STRING API] Jumlah baris interaksi mentah yang diterima: {len(raw_edges)}", log_fh)

        # Deduplikasi pasangan tak-berarah (undirected simple graph)
        unique_edges_0400 = {}
        for e in raw_edges:
            sid_a = e.get("stringId_A", "")
            sid_b = e.get("stringId_B", "")
            if sid_a not in string_to_gene or sid_b not in string_to_gene:
                continue
            g_a = string_to_gene[sid_a]
            g_b = string_to_gene[sid_b]
            if g_a == g_b:
                continue
            pair_key = tuple(sorted([g_a, g_b]))
            score = float(e.get("score", 0.0))
            if pair_key not in unique_edges_0400 or score > unique_edges_0400[pair_key]["combined_score"]:
                g1, g2 = pair_key
                unique_edges_0400[pair_key] = {
                    "protein1_stringId": gene_meta[g1]["stringId"],
                    "protein2_stringId": gene_meta[g2]["stringId"],
                    "gene1_id": g1,
                    "gene2_id": g2,
                    "symbol1": gene_meta[g1]["symbol"],
                    "symbol2": gene_meta[g2]["symbol"],
                    "direction1": gene_meta[g1]["direction"],
                    "direction2": gene_meta[g2]["direction"],
                    "combined_score": score,
                    "nscore": float(e.get("nscore", 0.0)),
                    "fscore": float(e.get("fscore", 0.0)),
                    "pscore": float(e.get("pscore", 0.0)),
                    "ascore": float(e.get("ascore", 0.0)),
                    "escore": float(e.get("escore", 0.0)),
                    "dscore": float(e.get("dscore", 0.0)),
                    "tscore": float(e.get("tscore", 0.0)),
                }

        edges_0400_list = sorted(unique_edges_0400.values(), key=lambda x: -x["combined_score"])
        edges_0700_list = [x for x in edges_0400_list if x["combined_score"] >= 0.700]

        log_print(f"[Tahap 3] Interaksi unik tak-berarah pada ambang >= 0.400 (Medium): {len(edges_0400_list)} edge", log_fh)
        log_print(f"[Tahap 3] Interaksi unik tak-berarah pada ambang >= 0.700 (High)  : {len(edges_0700_list)} edge", log_fh)

        # Simpan tabel interaksi
        edge_fields = [
            "protein1_stringId", "protein2_stringId", "gene1_id", "gene2_id",
            "symbol1", "symbol2", "direction1", "direction2",
            "combined_score", "nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"
        ]
        with open(os.path.join(RESULTS_DIR, "STRING_interactions_score0400_medium.csv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=edge_fields)
            writer.writeheader()
            writer.writerows(edges_0400_list)

        with open(os.path.join(RESULTS_DIR, "STRING_interactions_score0700_high.csv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=edge_fields)
            writer.writeheader()
            writer.writerows(edges_0700_list)

        # 3. Bangun graf NetworkX untuk 4 konfigurasi:
        # (A) Semua DEG (|log2FC|>=1.0) @ 0.400
        # (B) Semua DEG (|log2FC|>=1.0) @ 0.700
        # (C) Subset Ketat (|log2FC|>=1.5) @ 0.400
        # (D) Subset Ketat (|log2FC|>=1.5) @ 0.700
        def build_nx_graph(edge_list, min_abs_lfc=1.0):
            G = nx.Graph()
            for e in edge_list:
                g1, g2 = e["gene1_id"], e["gene2_id"]
                if gene_meta[g1]["abs_log2FC"] >= min_abs_lfc and gene_meta[g2]["abs_log2FC"] >= min_abs_lfc:
                    w = e["combined_score"]
                    G.add_edge(g1, g2, weight=w, distance=1.0 / w, **e)
            for u in G.nodes():
                G.nodes[u].update(gene_meta[u])
            return G

        G_0400_all = build_nx_graph(edges_0400_list, min_abs_lfc=1.0)
        G_0700_all = build_nx_graph(edges_0700_list, min_abs_lfc=1.0)
        G_0400_strict15 = build_nx_graph(edges_0400_list, min_abs_lfc=1.5)
        G_0700_strict15 = build_nx_graph(edges_0700_list, min_abs_lfc=1.5)

        n_deg_strict15 = sum(1 for r in all_deg_rows if r["abs_log2FC"] >= 1.5)
        n_mapped_strict15 = sum(1 for r in mapped_rows if r["abs_log2FC"] >= 1.5)

        stats_0400_all, GCC_0400_all = compute_network_stats(
            G_0400_all, len(all_deg_rows), len(mapped_rows),
            "1. Semua DEG (|log2FC|>=1.0) - Skor >= 0.400 (Medium, Default Utama)"
        )
        stats_0700_all, GCC_0700_all = compute_network_stats(
            G_0700_all, len(all_deg_rows), len(mapped_rows),
            "2. Semua DEG (|log2FC|>=1.0) - Skor >= 0.700 (High, Sensitivitas Skor)"
        )
        stats_0400_s15, GCC_0400_s15 = compute_network_stats(
            G_0400_strict15, n_deg_strict15, n_mapped_strict15,
            "3. Subset DEG Ketat (|log2FC|>=1.5) - Skor >= 0.400 (Sensitivitas Kepadatan)"
        )
        stats_0700_s15, GCC_0700_s15 = compute_network_stats(
            G_0700_strict15, n_deg_strict15, n_mapped_strict15,
            "4. Subset DEG Ketat (|log2FC|>=1.5) - Skor >= 0.700 (Ketat & High Confidence)"
        )

        all_stats = [stats_0400_all, stats_0700_all, stats_0400_s15, stats_0700_s15]
        stats_csv = os.path.join(RESULTS_DIR, "STRING_network_topology_statistics.csv")
        with open(stats_csv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(stats_0400_all.keys()))
            writer.writeheader()
            writer.writerows(all_stats)

        log_print("\n[Tahap 3] Ringkasan Statistik Topologi Jaringan STRING:", log_fh)
        for st in all_stats:
            log_print(
                f"   * {st['Konfigurasi_Jaringan']}:\n"
                f"     - Protein Terpetakan: {st['Protein_Terpetakan_STRING']} | Node Aktif: {st['Node_Aktif_Berinteraksi']} | Node Terisolasi: {st['Node_Terisolasi_Degree0']}\n"
                f"     - Jumlah Edge: {st['Jumlah_Edge']} | Avg Degree (Semua Terpetakan): {st['Rata2_Degree_Semua_Protein_Terpetakan']} | Avg Degree (Node Aktif): {st['Rata2_Degree_Node_Aktif']}\n"
                f"     - Komponen Terhubung (Aktif): {st['Jumlah_Komponen_Terhubung_Aktif']} | GCC: {st['Node_GCC_Komponen_Terbesar']} node ({st['Persen_Node_Aktif_di_GCC']}%), {st['Edge_GCC_Komponen_Terbesar']} edge, Clustering={st['Clustering_Coefficient_GCC']}, ASPL={st['Rata2_Jarak_Terpendek_GCC']}",
                log_fh
            )

        # 4. Tahap 4: Deteksi Komunitas Louvain pada GCC Utama (0.400)
        log_print("\n[Tahap 4] Menjalankan deteksi komunitas Louvain pada GCC Utama (Skor >= 0.400)...", log_fh)
        louvain_comms = nx.community.louvain_communities(GCC_0400_all, weight="weight", seed=SEED)
        louvain_comms = sorted(louvain_comms, key=len, reverse=True)
        modularity_q = nx.community.modularity(GCC_0400_all, louvain_comms, weight="weight")
        log_print(f"[Louvain] Jumlah komunitas/modul pada GCC (547+ node): {len(louvain_comms)} modul | Modularitas Q = {modularity_q:.4f}", log_fh)

        node_to_module = {}
        for m_idx, comm in enumerate(louvain_comms, start=1):
            mod_name = f"M{m_idx}"
            for u in comm:
                node_to_module[u] = mod_name

        # 5. Tahap 4 & 5: Perhitungan Centrality & Composite Score pada 5 Skema untuk Uji Stabilitas
        log_print("\n[Tahap 4 & 5] Menghitung Degree, Betweenness, Closeness, Eigenvector Centrality & Composite Score...", log_fh)
        cent_gcc_0400_unw = compute_centralities_and_composite(GCC_0400_all, weighted=False)
        cent_gcc_0400_w   = compute_centralities_and_composite(GCC_0400_all, weighted=True)
        cent_full_0400_unw = compute_centralities_and_composite(G_0400_all, weighted=False)
        cent_gcc_0700_unw = compute_centralities_and_composite(GCC_0700_all, weighted=False)
        cent_gcc_0700_w   = compute_centralities_and_composite(GCC_0700_all, weighted=True)

        # Gabungkan informasi lengkap untuk seluruh node di GCC 0.400 (Utama)
        gcc_ranked_rows = []
        for gid, cinfo in sorted(cent_gcc_0400_unw.items(), key=lambda x: x[1]["rank"]):
            meta = gene_meta[gid]
            w_info = cent_gcc_0400_w.get(gid, {})
            f_info = cent_full_0400_unw.get(gid, {})
            c7_unw = cent_gcc_0700_unw.get(gid, {})
            c7_w   = cent_gcc_0700_w.get(gid, {})

            row = {
                "rank_primary": cinfo["rank"],
                "gene_id": gid,
                "symbol": meta["symbol"],
                "stringId": meta["stringId"],
                "direction": meta["direction"],
                "regulation_label": "Naik pada Kitin (Up)" if meta["direction"] == "Up" else "Turun pada Kitin (Down)",
                "log2FC": round(meta["log2FC"], 4),
                "padj": meta["padj"],
                "baseMean": round(meta["baseMean"], 2),
                "louvain_module": node_to_module.get(gid, "Non-GCC"),
                "degree_raw": cinfo["degree_raw"],
                "degree_centrality": round(cinfo["degree_centrality"], 6),
                "betweenness_centrality": round(cinfo["betweenness_centrality"], 6),
                "closeness_centrality": round(cinfo["closeness_centrality"], 6),
                "eigenvector_centrality": round(cinfo["eigenvector_centrality"], 6),
                "z_degree": round(cinfo["z_degree"], 4),
                "z_betweenness": round(cinfo["z_betweenness"], 4),
                "z_closeness": round(cinfo["z_closeness"], 4),
                "z_eigenvector": round(cinfo["z_eigenvector"], 4),
                "composite_score": round(cinfo["composite_score"], 4),
                "rank_0400_weighted": w_info.get("rank", "NA"),
                "composite_0400_weighted": round(w_info["composite_score"], 4) if "composite_score" in w_info else "NA",
                "rank_0400_full_active": f_info.get("rank", "NA"),
                "rank_0700_unweighted": c7_unw.get("rank", "Non-GCC_0700"),
                "degree_raw_0700": c7_unw.get("degree_raw", 0),
                "composite_0700_unweighted": round(c7_unw["composite_score"], 4) if "composite_score" in c7_unw else "NA",
                "rank_0700_weighted": c7_w.get("rank", "Non-GCC_0700"),
                "description": meta["description"],
                "string_annotation": meta["string_annotation"]
            }
            gcc_ranked_rows.append(row)

        top10_rows = gcc_ranked_rows[:10]

        log_print("\n[Tahap 5] TOP 10 KANDIDAT GEN PRIORITAS (GCC 0.400 Unweighted, Z-Score STDEV.S):", log_fh)
        for r in top10_rows:
            log_print(
                f"   #{r['rank_primary']:2d} | {r['gene_id']} ({r['symbol']:10s}) | Dir={r['direction']:4s} | log2FC={r['log2FC']:+6.3f} | padj={r['padj']:.2e} | "
                f"Deg={r['degree_raw']:3d} (C_D={r['degree_centrality']:.4f}, C_B={r['betweenness_centrality']:.4f}, C_C={r['closeness_centrality']:.4f}, C_E={r['eigenvector_centrality']:.4f}) | "
                f"Composite={r['composite_score']:+6.3f} | Mod={r['louvain_module']} | Ranks[0.4W={r['rank_0400_weighted']}, Full={r['rank_0400_full_active']}, 0.7U={r['rank_0700_unweighted']}, 0.7W={r['rank_0700_weighted']}]",
                log_fh
            )

        # Simpan tabel centrality lengkap
        cent_fields = list(gcc_ranked_rows[0].keys())
        with open(os.path.join(RESULTS_DIR, "PPI_centrality_GCC_score0400_complete.csv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cent_fields)
            writer.writeheader()
            writer.writerows(gcc_ranked_rows)

        with open(os.path.join(RESULTS_DIR, "top_candidates.csv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cent_fields)
            writer.writeheader()
            writer.writerows(top10_rows)

        # Simpan tabel ringkasan komunitas Louvain
        comm_rows = []
        for m_idx, comm in enumerate(louvain_comms, start=1):
            mod_name = f"M{m_idx}"
            sub_nodes = sorted(list(comm), key=lambda u: cent_gcc_0400_unw[u]["composite_score"], reverse=True)
            n_up = sum(1 for u in sub_nodes if gene_meta[u]["direction"] == "Up")
            n_down = len(sub_nodes) - n_up
            top_hubs = [f"{gene_meta[u]['symbol']} ({u}, k={cent_gcc_0400_unw[u]['degree_raw']})" for u in sub_nodes[:5]]
            comm_rows.append({
                "Module_ID": mod_name,
                "Jumlah_Node": len(sub_nodes),
                "Node_Naik_Up": n_up,
                "Node_Turun_Down": n_down,
                "Persen_Up": round(n_up / len(sub_nodes) * 100.0, 1),
                "Hub_Utama_Modul": "; ".join(top_hubs)
            })

        with open(os.path.join(RESULTS_DIR, "PPI_Louvain_communities_summary.csv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(comm_rows[0].keys()))
            writer.writeheader()
            writer.writerows(comm_rows)

        # 6. Uji Kestabilan Peringkat (Rank Stability Analysis)
        schemes = {
            "0.400_GCC_Unweighted (Utama)": cent_gcc_0400_unw,
            "0.400_GCC_Weighted": cent_gcc_0400_w,
            "0.400_FullActive_Unweighted": cent_full_0400_unw,
            "0.700_GCC_Unweighted": cent_gcc_0700_unw,
            "0.700_GCC_Weighted": cent_gcc_0700_w
        }
        top10_sets = {}
        top10_ordered = {}
        for s_name, s_dict in schemes.items():
            ordered = [u for u, _ in sorted(s_dict.items(), key=lambda x: x[1]["rank"])[:10]]
            top10_ordered[s_name] = ordered
            top10_sets[s_name] = set(ordered)

        ref_name = "0.400_GCC_Unweighted (Utama)"
        ref_set = top10_sets[ref_name]
        ref_dict = schemes[ref_name]

        stability_summary = []
        for s_name, s_dict in schemes.items():
            common_nodes = sorted(list(set(ref_dict.keys()).intersection(set(s_dict.keys()))))
            r1 = [ref_dict[u]["composite_score"] for u in common_nodes]
            r2 = [s_dict[u]["composite_score"] for u in common_nodes]
            rho, p_rho = spearmanr(r1, r2)
            tau, p_tau = kendalltau(r1, r2)
            overlap = ref_set.intersection(top10_sets[s_name])
            jaccard = len(overlap) / len(ref_set.union(top10_sets[s_name]))
            top10_syms = [f"#{i+1} {gene_meta[u]['symbol']} ({u})" for i, u in enumerate(top10_ordered[s_name])]
            overlap_syms = [gene_meta[u]["symbol"] for u in top10_ordered[ref_name] if u in overlap]
            stability_summary.append({
                "Skema_Pembanding": s_name,
                "Jumlah_Node_Bersama": len(common_nodes),
                "Irisan_Top10_dengan_Utama": f"{len(overlap)}/10",
                "Jaccard_Index_Top10": round(jaccard, 4),
                "Spearman_Rho_Composite": round(float(rho), 4),
                "Kendall_Tau_Composite": round(float(tau), 4),
                "Gen_Top10_Konsisten": ", ".join(overlap_syms),
                "Urutan_Top10_Pada_Skema_Ini": "; ".join(top10_syms)
            })

        with open(os.path.join(RESULTS_DIR, "PPI_top10_stability_comparison.csv"), "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(stability_summary[0].keys()))
            writer.writeheader()
            writer.writerows(stability_summary)

        log_print("\n[Tahap 5] Hasil Uji Kestabilan Peringkat terhadap Skema Utama (0.400 GCC Unweighted):", log_fh)
        for st in stability_summary:
            log_print(
                f"   * {st['Skema_Pembanding']:30s} | Irisan Top10={st['Irisan_Top10_dengan_Utama']} (Jaccard={st['Jaccard_Index_Top10']:.3f}) | "
                f"Spearman rho={st['Spearman_Rho_Composite']:.4f} | Kendall tau={st['Kendall_Tau_Composite']:.4f}\n"
                f"     Top 10: {st['Urutan_Top10_Pada_Skema_Ini']}",
                log_fh
            )

        # 7. Simpan File Excel `results/top_candidates.xlsx`
        xlsx_path = os.path.join(RESULTS_DIR, "top_candidates.xlsx")
        wb = openpyxl.Workbook()
        # Sheet 1: Top 10 Candidates
        ws1 = wb.active
        ws1.title = "Top10_Candidates"

        header_font = Font(bold=True, color="FFFFFF", name="Calibri", size=11)
        header_fill = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
        up_fill = PatternFill(start_color="FCE4D6", end_color="FCE4D6", fill_type="solid")
        down_fill = PatternFill(start_color="DDEBF7", end_color="DDEBF7", fill_type="solid")
        thin_border = Border(
            left=Side(style="thin", color="D9D9D9"),
            right=Side(style="thin", color="D9D9D9"),
            top=Side(style="thin", color="D9D9D9"),
            bottom=Side(style="thin", color="D9D9D9")
        )

        def write_sheet(ws, rows, fieldnames):
            ws.append(fieldnames)
            for col_idx in range(1, len(fieldnames) + 1):
                cell = ws.cell(row=1, column=col_idx)
                cell.font = header_font
                cell.fill = header_fill
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            for r_idx, r in enumerate(rows, start=2):
                ws.append([r.get(k, "") for k in fieldnames])
                fill_to_use = up_fill if r.get("direction") == "Up" else (down_fill if r.get("direction") == "Down" else None)
                for col_idx in range(1, len(fieldnames) + 1):
                    c = ws.cell(row=r_idx, column=col_idx)
                    c.border = thin_border
                    if fill_to_use and col_idx <= 8:
                        c.fill = fill_to_use
            for col_idx, col_name in enumerate(fieldnames, start=1):
                col_letter = get_column_letter(col_idx)
                max_len = max(len(str(col_name)), *(len(str(r.get(col_name, ""))) for r in rows[:50]))
                ws.column_dimensions[col_letter].width = min(max(max_len + 3, 12), 45)
            ws.freeze_panes = "A2"

        write_sheet(ws1, top10_rows, cent_fields)

        ws2 = wb.create_sheet(title="All_GCC_Nodes_Ranked")
        write_sheet(ws2, gcc_ranked_rows, cent_fields)

        ws3 = wb.create_sheet(title="Rank_Stability_Comparison")
        write_sheet(ws3, stability_summary, list(stability_summary[0].keys()))

        ws4 = wb.create_sheet(title="Network_Topology_Stats")
        write_sheet(ws4, all_stats, list(stats_0400_all.keys()))

        wb.save(xlsx_path)
        log_print(f"\n[Output Excel] Berkas 'results/top_candidates.xlsx' berhasil disimpan ({xlsx_path}).", log_fh)

        # =========================================================================
        # 8. VISUALISASI JARINGAN PPI & CENTRALITY (PNG 300 DPI + PDF + SVG)
        # =========================================================================
        log_print("\n[Visualisasi] Membuat grafik jaringan PPI, komunitas, distribusi degree, dan profil centrality Top 10...", log_fh)
        top10_ids = set(r["gene_id"] for r in top10_rows)

        # --- GAMBAR PPI 03: Jaringan PPI Utama (Panel A: GCC 0.400 & Panel B: GCC 0.700 High Confidence) ---
        fig, axes = plt.subplots(1, 2, figsize=(18, 8.5), dpi=300)
        fig.suptitle(
            "Jaringan Interaksi Protein-Protein (PPI) STRING DEG Respons Kitin (C6 1h vs Kontrol H2O 1h)\n"
            "Organisme: Arabidopsis thaliana (NCBI TaxID: 3702) | Warna Node: Arah Regulasi Kontras Perlakuan vs Kontrol",
            fontsize=13, fontweight="bold", y=0.98
        )

        # Panel A: GCC 0.400 (Default Medium Confidence)
        ax = axes[0]
        pos_0400 = nx.spring_layout(GCC_0400_all, k=0.26, iterations=80, seed=SEED, weight="weight")
        # Gambar edge tipis
        nx.draw_networkx_edges(
            GCC_0400_all, pos_0400, ax=ax,
            edge_color="#B0B0B0", alpha=0.22, width=0.5
        )
        # Gambar node non-Top10
        non_top_0400 = [u for u in GCC_0400_all.nodes() if u not in top10_ids]
        node_colors_non_top = [COLOR_UP if gene_meta[u]["direction"] == "Up" else COLOR_DOWN for u in non_top_0400]
        node_sizes_non_top = [18 + 3.2 * cent_gcc_0400_unw[u]["degree_raw"] for u in non_top_0400]
        nx.draw_networkx_nodes(
            GCC_0400_all, pos_0400, nodelist=non_top_0400, ax=ax,
            node_color=node_colors_non_top, node_size=node_sizes_non_top,
            alpha=0.82, linewidths=0.3, edgecolors="white"
        )
        # Gambar node Top 10 dengan border emas/hitam tegas
        top_list_0400 = [r["gene_id"] for r in top10_rows]
        top_colors_0400 = [COLOR_UP if gene_meta[u]["direction"] == "Up" else COLOR_DOWN for u in top_list_0400]
        top_sizes_0400 = [160 + 5.5 * cent_gcc_0400_unw[u]["degree_raw"] for u in top_list_0400]
        nx.draw_networkx_nodes(
            GCC_0400_all, pos_0400, nodelist=top_list_0400, ax=ax,
            node_color=top_colors_0400, node_size=top_sizes_0400,
            alpha=1.0, linewidths=1.8, edgecolors="#111111"
        )
        # Label Top 10
        for r in top10_rows:
            u = r["gene_id"]
            x, y = pos_0400[u]
            ax.annotate(
                f"#{r['rank_primary']} {r['symbol']}",
                xy=(x, y), xytext=(x + 0.03, y + 0.03),
                fontsize=8.5, fontweight="bold", color="#111111",
                bbox=dict(boxstyle="round,pad=0.18", fc="#FFFDE7", ec="#333333", alpha=0.9, lw=0.7)
            )
        ax.set_title(
            f"A. Komponen Terhubung Terbesar (GCC) Skor >= 0.400 (Medium Confidence)\n"
            f"({GCC_0400_all.number_of_nodes()} Node, {GCC_0400_all.number_of_edges()} Edge, Avg Degree = {2*GCC_0400_all.number_of_edges()/GCC_0400_all.number_of_nodes():.2f})",
            fontsize=11, fontweight="bold"
        )
        ax.axis("off")

        # Panel B: GCC 0.700 (High Confidence)
        ax2 = axes[1]
        pos_0700 = nx.spring_layout(GCC_0700_all, k=0.34, iterations=90, seed=SEED, weight="weight")
        nx.draw_networkx_edges(
            GCC_0700_all, pos_0700, ax=ax2,
            edge_color="#888888", alpha=0.38, width=0.8
        )
        non_top_0700 = [u for u in GCC_0700_all.nodes() if u not in top10_ids]
        node_colors_0700 = [COLOR_UP if gene_meta[u]["direction"] == "Up" else COLOR_DOWN for u in non_top_0700]
        node_sizes_0700 = [25 + 6.0 * cent_gcc_0700_unw[u]["degree_raw"] for u in non_top_0700]
        nx.draw_networkx_nodes(
            GCC_0700_all, pos_0700, nodelist=non_top_0700, ax=ax2,
            node_color=node_colors_0700, node_size=node_sizes_0700,
            alpha=0.85, linewidths=0.4, edgecolors="white"
        )
        top_in_0700 = [u for u in top_list_0400 if u in GCC_0700_all]
        top_colors_0700 = [COLOR_UP if gene_meta[u]["direction"] == "Up" else COLOR_DOWN for u in top_in_0700]
        top_sizes_0700 = [170 + 8.0 * cent_gcc_0700_unw[u]["degree_raw"] for u in top_in_0700]
        nx.draw_networkx_nodes(
            GCC_0700_all, pos_0700, nodelist=top_in_0700, ax=ax2,
            node_color=top_colors_0700, node_size=top_sizes_0700,
            alpha=1.0, linewidths=1.8, edgecolors="#111111"
        )
        for r in top10_rows:
            u = r["gene_id"]
            if u in pos_0700:
                x, y = pos_0700[u]
                ax2.annotate(
                    f"#{r['rank_primary']} {r['symbol']}",
                    xy=(x, y), xytext=(x + 0.03, y + 0.03),
                    fontsize=8.5, fontweight="bold", color="#111111",
                    bbox=dict(boxstyle="round,pad=0.18", fc="#FFFDE7", ec="#333333", alpha=0.9, lw=0.7)
                )
        ax2.set_title(
            f"B. Komponen Terhubung Terbesar (GCC) Skor >= 0.700 (High Confidence)\n"
            f"({GCC_0700_all.number_of_nodes()} Node, {GCC_0700_all.number_of_edges()} Edge, Avg Degree = {2*GCC_0700_all.number_of_edges()/max(1,GCC_0700_all.number_of_nodes()):.2f})",
            fontsize=11, fontweight="bold"
        )
        ax2.axis("off")

        legend_elems = [
            Line2D([0], [0], marker="o", color="w", label="DEG Naik pada Kitin C6 (log2FC >= +1.0)", markerfacecolor=COLOR_UP, markersize=10),
            Line2D([0], [0], marker="o", color="w", label="DEG Turun pada Kitin C6 (log2FC <= -1.0)", markerfacecolor=COLOR_DOWN, markersize=10),
            Line2D([0], [0], marker="o", color="w", label="Top 10 Kandidat Hub Prioritas (Border Hitam)", markerfacecolor=COLOR_UP, markeredgecolor="#111111", markeredgewidth=1.8, markersize=12),
        ]
        fig.legend(handles=legend_elems, loc="lower center", ncol=3, frameon=True, fontsize=10, bbox_to_anchor=(0.5, 0.01))
        plt.tight_layout(rect=[0, 0.05, 1, 0.93])
        save_fig_multi(fig, "Fig_PPI_03_STRING_PPI_Network_Overview")
        plt.close(fig)

        # --- GAMBAR PPI 04: Komunitas Louvain & Sub-Jaringan Inti Top 10 Hub + Tetangga Langsung ---
        fig, axes = plt.subplots(1, 2, figsize=(18, 8.5), dpi=300)
        fig.suptitle(
            "Modularitas Komunitas Louvain & Sub-Jaringan Interaksi Inti di Sekitar Top 10 Kandidat Hub\n"
            "Kontras: Perlakuan Kitin C6 1 jam vs Kontrol H2O 1 jam (Arabidopsis thaliana, GSE133053)",
            fontsize=13, fontweight="bold", y=0.98
        )

        # Panel A: Pewarnaan 6 Modul Louvain Terbesar pada GCC 0.400
        ax = axes[0]
        palette_mod = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#999999"]
        mod_colors = []
        for u in GCC_0400_all.nodes():
            m_str = node_to_module.get(u, "M99")
            m_num = int(m_str[1:]) if m_str[1:].isdigit() else 99
            mod_colors.append(palette_mod[m_num - 1] if m_num <= 6 else "#BBBBBB")

        nx.draw_networkx_edges(GCC_0400_all, pos_0400, ax=ax, edge_color="#CCCCCC", alpha=0.2, width=0.5)
        nx.draw_networkx_nodes(
            GCC_0400_all, pos_0400, ax=ax,
            node_color=mod_colors,
            node_size=[22 + 3.5 * cent_gcc_0400_unw[u]["degree_raw"] for u in GCC_0400_all.nodes()],
            alpha=0.88, linewidths=0.3, edgecolors="white"
        )
        ax.set_title(
            f"A. Deteksi Komunitas Louvain pada GCC (Skor >= 0.400)\n"
            f"Total {len(louvain_comms)} Komunitas | Modularitas Q = {modularity_q:.3f}",
            fontsize=11, fontweight="bold"
        )
        ax.axis("off")
        mod_patches = [
            mpatches.Patch(color=palette_mod[i], label=f"Modul M{i+1} (n={len(louvain_comms[i])} gen, Top: {gene_meta[sorted(list(louvain_comms[i]), key=lambda x: cent_gcc_0400_unw[x]['composite_score'], reverse=True)[0]]['symbol']})")
            for i in range(min(6, len(louvain_comms)))
        ]
        ax.legend(handles=mod_patches, loc="lower left", fontsize=8.2, frameon=True, title="6 Modul Louvain Terbesar")

        # Panel B: Sub-jaringan induksi di antara Top 35 node berdasarkan Composite Score (agar sangat terbaca)
        ax2 = axes[1]
        top35_nodes = [r["gene_id"] for r in gcc_ranked_rows[:35]]
        G_sub35 = GCC_0400_all.subgraph(top35_nodes).copy()
        pos_sub35 = nx.kamada_kawai_layout(G_sub35)
        edge_widths_sub = [1.2 + 2.5 * (d.get("weight", 0.4) - 0.4) for _, _, d in G_sub35.edges(data=True)]
        nx.draw_networkx_edges(G_sub35, pos_sub35, ax=ax2, edge_color="#777777", alpha=0.45, width=edge_widths_sub)

        sub_colors = [COLOR_UP if gene_meta[u]["direction"] == "Up" else COLOR_DOWN for u in G_sub35.nodes()]
        sub_sizes = [260 if u in top10_ids else 110 for u in G_sub35.nodes()]
        sub_borders = ["#111111" if u in top10_ids else "white" for u in G_sub35.nodes()]
        sub_lw = [2.0 if u in top10_ids else 0.7 for u in G_sub35.nodes()]
        nx.draw_networkx_nodes(
            G_sub35, pos_sub35, ax=ax2,
            node_color=sub_colors, node_size=sub_sizes,
            edgecolors=sub_borders, linewidths=sub_lw, alpha=0.95
        )
        for u in G_sub35.nodes():
            x, y = pos_sub35[u]
            sym = gene_meta[u]["symbol"]
            is_t10 = u in top10_ids
            rk = cent_gcc_0400_unw[u]["rank"]
            lbl = f"#{rk} {sym}" if is_t10 else sym
            ax2.text(
                x, y + 0.045, lbl,
                fontsize=8.5 if is_t10 else 7.2,
                fontweight="bold" if is_t10 else "normal",
                ha="center", va="bottom",
                color="#000000" if is_t10 else "#333333",
                bbox=dict(boxstyle="round,pad=0.12", fc="#FFF9C4" if is_t10 else "white", ec="none", alpha=0.82)
            )
        ax2.set_title(
            f"B. Sub-Jaringan Inti Top 35 Hub Prioritas ({G_sub35.number_of_nodes()} Node, {G_sub35.number_of_edges()} Interaksi Langsung)\n"
            f"Menampilkan Keterhubungan Erat Antar-Kandidat Top 10 (Kotak Kuning)",
            fontsize=11, fontweight="bold"
        )
        ax2.axis("off")

        plt.tight_layout(rect=[0, 0.02, 1, 0.93])
        save_fig_multi(fig, "Fig_PPI_04_Louvain_Community_and_Hub_Subnetwork")
        plt.close(fig)

        # --- GAMBAR PPI 05: Distribusi Degree & Karakteristik Topologi Scale-Free (0.400 vs 0.700) ---
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), dpi=300)
        fig.suptitle(
            "Distribusi Degree, Uji Topologi Scale-Free (Log-Log), dan Sensitivitas Kepadatan Jaringan PPI STRING\n"
            "Perlakuan Kitin C6 1 jam vs Kontrol H2O 1 jam (GSE133053)",
            fontsize=12.5, fontweight="bold", y=1.01
        )

        # Panel A: Histogram Degree pada GCC 0.400
        ax = axes[0]
        degs_0400 = [cent_gcc_0400_unw[u]["degree_raw"] for u in GCC_0400_all.nodes()]
        ax.hist(degs_0400, bins=25, color="#0072B2", edgecolor="white", alpha=0.85)
        ax.axvline(np.mean(degs_0400), color="#D55E00", linestyle="--", linewidth=2, label=f"Mean Degree = {np.mean(degs_0400):.2f}")
        ax.axvline(np.median(degs_0400), color="#009E73", linestyle=":", linewidth=2, label=f"Median Degree = {np.median(degs_0400):.1f}")
        # Tandai Top 3 hub
        for r in top10_rows[:3]:
            ax.annotate(
                f"#{r['rank_primary']} {r['symbol']}\n(k={r['degree_raw']})",
                xy=(r["degree_raw"], 3), xytext=(r["degree_raw"] - 8, 28 + 12 * (r["rank_primary"] % 2)),
                arrowprops=dict(arrowstyle="->", color="#111111", lw=1.0),
                fontsize=8.5, fontweight="bold"
            )
        ax.set_xlabel("Degree (Jumlah Interaksi per Protein, k)", fontsize=10, fontweight="bold")
        ax.set_ylabel("Jumlah Protein (Frekuensi)", fontsize=10, fontweight="bold")
        ax.set_title("A. Distribusi Degree pada GCC (Skor >= 0.400)", fontsize=10.5, fontweight="bold")
        ax.legend(frameon=True, fontsize=9)
        ax.grid(axis="y", alpha=0.25)

        # Panel B: Plot Log-Log P(k) vs k (Hukum Pangkat / Scale-free) pada 0.400 vs 0.700
        ax2 = axes[1]
        for deg_list, lbl, col, mk in [
            (degs_0400, "Skor >= 0.400 (Medium)", "#0072B2", "o"),
            ([cent_gcc_0700_unw[u]["degree_raw"] for u in GCC_0700_all.nodes()], "Skor >= 0.700 (High)", "#D55E00", "s")
        ]:
            vals, counts = np.unique(deg_list, return_counts=True)
            probs = counts / np.sum(counts)
            mask = vals > 0
            log_k = np.log10(vals[mask])
            log_p = np.log10(probs[mask])
            slope, intercept = np.polyfit(log_k, log_p, 1)
            r_corr = np.corrcoef(log_k, log_p)[0, 1]
            ax2.scatter(vals[mask], probs[mask], color=col, marker=mk, s=40, alpha=0.8, label=f"{lbl} (gamma={-slope:.2f}, RÂ²={r_corr**2:.2f})")
            k_line = np.linspace(min(vals[mask]), max(vals[mask]), 50)
            ax2.plot(k_line, (10**intercept) * (k_line**slope), color=col, linestyle="--", linewidth=1.5)

        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlabel("Degree (k) [Skala Log10]", fontsize=10, fontweight="bold")
        ax2.set_ylabel("Probabilitas P(k) [Skala Log10]", fontsize=10, fontweight="bold")
        ax2.set_title("B. Fit Hukum Pangkat P(k) ~ k^(-gamma)", fontsize=10.5, fontweight="bold")
        ax2.legend(frameon=True, fontsize=8.5)
        ax2.grid(True, which="both", alpha=0.25)

        # Panel C: Perbandingan Jumlah Node Aktif, GCC, dan Isolat antar-4 Konfigurasi
        ax3 = axes[2]
        cfg_labels = ["All DEG\n(>=0.400)", "All DEG\n(>=0.700)", "|log2FC|>=1.5\n(>=0.400)", "|log2FC|>=1.5\n(>=0.700)"]
        gcc_vals = [st["Node_GCC_Komponen_Terbesar"] for st in all_stats]
        other_active_vals = [st["Node_Aktif_Berinteraksi"] - st["Node_GCC_Komponen_Terbesar"] for st in all_stats]
        isol_vals = [st["Node_Terisolasi_Degree0"] for st in all_stats]

        x_pos = np.arange(len(cfg_labels))
        b1 = ax3.bar(x_pos, gcc_vals, width=0.55, color="#0072B2", label="Node di Komponen Terbesar (GCC)")
        b2 = ax3.bar(x_pos, other_active_vals, width=0.55, bottom=gcc_vals, color="#56B4E9", label="Node di Komponen Kecil")
        b3 = ax3.bar(x_pos, isol_vals, width=0.55, bottom=np.array(gcc_vals) + np.array(other_active_vals), color="#CCCCCC", label="Node Terisolasi (Degree = 0)")

        for idx_b, st in enumerate(all_stats):
            ax3.text(
                idx_b, st["Node_GCC_Komponen_Terbesar"] / 2,
                f"GCC: {st['Node_GCC_Komponen_Terbesar']}\n({st['Jumlah_Edge']} edge)",
                ha="center", va="center", color="white", fontsize=8, fontweight="bold"
            )
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(cfg_labels, fontsize=9, fontweight="bold")
        ax3.set_ylabel("Jumlah Protein Terpetakan STRING", fontsize=10, fontweight="bold")
        ax3.set_title("C. Komposisi Komponen Jaringan per Ambang", fontsize=10.5, fontweight="bold")
        ax3.legend(frameon=True, fontsize=8.5, loc="upper right")
        ax3.grid(axis="y", alpha=0.25)

        plt.tight_layout()
        save_fig_multi(fig, "Fig_PPI_05_Degree_Distribution_and_Topology")
        plt.close(fig)

        # --- GAMBAR PPI 06: Skor Prioritas Komposit Top 10 & Uji Stabilitas Peringkat ---
        fig, axes = plt.subplots(1, 2, figsize=(17, 7.2), dpi=300, gridspec_kw={"width_ratios": [1.25, 1.0]})
        fig.suptitle(
            "Prioritisasi Top 10 Kandidat Gen Hub Berdasarkan Skor Komposit Centrality & Uji Kestabilan Peringkat\n"
            "Kontras: Perlakuan Kitin C6 1 jam vs Kontrol H2O 1 jam (GSE133053)",
            fontsize=12.5, fontweight="bold", y=0.99
        )

        # Panel A: Horizontal Bar / Dot Plot Composite Score & Kontribusi 4 Z-Score
        ax = axes[0]
        y_pos = np.arange(len(top10_rows))[::-1]
        labels_t10 = [f"#{r['rank_primary']} {r['symbol']} ({r['gene_id']})\n[log2FC={r['log2FC']:+.2f}, k={r['degree_raw']}]" for r in top10_rows]
        comp_vals = [r["composite_score"] for r in top10_rows]
        bar_cols = [COLOR_UP if r["direction"] == "Up" else COLOR_DOWN for r in top10_rows]

        bars = ax.barh(y_pos, comp_vals, height=0.62, color=bar_cols, alpha=0.85, edgecolor="#222222", linewidth=0.9, label="Composite Score (Rata-rata 4 Z-score)")

        # Plot titik masing-masing Z-score
        for z_key, z_lbl, z_mk, z_col in [
            ("z_degree", "Z-Degree", "o", "#000000"),
            ("z_betweenness", "Z-Betweenness", "s", "#009E73"),
            ("z_closeness", "Z-Closeness", "^", "#56B4E9"),
            ("z_eigenvector", "Z-Eigenvector", "D", "#F0E442")
        ]:
            ax.scatter(
                [r[z_key] for r in top10_rows], y_pos,
                marker=z_mk, s=52, color=z_col, edgecolors="#111111", linewidths=0.8, zorder=5, label=z_lbl
            )

        for idx_r, r in enumerate(top10_rows):
            ax.text(
                r["composite_score"] + 0.12, y_pos[idx_r],
                f"{r['composite_score']:.2f}",
                va="center", ha="left", fontsize=9, fontweight="bold", color="#111111"
            )

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels_t10, fontsize=9, fontweight="bold")
        ax.set_xlabel("Nilai Z-Score (STDEV.S) / Skor Prioritas Komposit", fontsize=10, fontweight="bold")
        ax.set_title("A. Profil 4 Z-Score Centrality & Skor Komposit Top 10 Hub", fontsize=11, fontweight="bold")
        ax.legend(loc="lower right", frameon=True, fontsize=8.5, ncol=2)
        ax.grid(axis="x", alpha=0.3)
        ax.set_xlim(0, max(max(r["z_betweenness"] for r in top10_rows), max(comp_vals)) + 1.2)

        # Panel B: Heatmap Peringkat Top 10 melintasi 5 Skema Uji Stabilitas
        ax2 = axes[1]
        scheme_short = [
            "0.400 GCC\nUnweighted\n(Utama)",
            "0.400 GCC\nWeighted\n(Skor)",
            "0.400 Full\nActive\nUnweighted",
            "0.700 GCC\nUnweighted\n(Ketat)",
            "0.700 GCC\nWeighted\n(Ketat)"
        ]
        rank_matrix = np.zeros((10, 5), dtype=float)
        for i_r, r in enumerate(top10_rows):
            gid = r["gene_id"]
            rank_matrix[i_r, 0] = cent_gcc_0400_unw[gid]["rank"]
            rank_matrix[i_r, 1] = cent_gcc_0400_w[gid]["rank"]
            rank_matrix[i_r, 2] = cent_full_0400_unw[gid]["rank"]
            rank_matrix[i_r, 3] = cent_gcc_0700_unw[gid]["rank"] if gid in cent_gcc_0700_unw else np.nan
            rank_matrix[i_r, 4] = cent_gcc_0700_w[gid]["rank"] if gid in cent_gcc_0700_w else np.nan

        # Gunakan skala warna terbalik (peringkat 1 paling gelap/kuat)
        masked_mat = np.ma.masked_invalid(rank_matrix)
        im = ax2.imshow(np.clip(masked_mat, 1, 30), cmap="YlGnBu_r", aspect="auto", vmin=1, vmax=25)
        for i_r in range(10):
            for j_c in range(5):
                val = rank_matrix[i_r, j_c]
                txt = f"#{int(val)}" if not np.isnan(val) else "Isolat"
                col_txt = "white" if (not np.isnan(val) and val <= 10) else "#111111"
                ax2.text(j_c, i_r, txt, ha="center", va="center", color=col_txt, fontsize=9.5, fontweight="bold")

        ax2.set_xticks(np.arange(5))
        ax2.set_xticklabels(scheme_short, fontsize=8.8, fontweight="bold")
        ax2.set_yticks(np.arange(10))
        ax2.set_yticklabels([f"#{r['rank_primary']} {r['symbol']}" for r in top10_rows], fontsize=9.5, fontweight="bold")
        ax2.set_title("B. Stabilitas Peringkat Top 10 pada Berbagai Skema", fontsize=11, fontweight="bold")
        cbar = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
        cbar.set_label("Peringkat Komposit (1 = Tertinggi)", fontsize=9, fontweight="bold")

        plt.tight_layout(rect=[0, 0.01, 1, 0.94])
        save_fig_multi(fig, "Fig_PPI_06_Top10_Centrality_and_Composite_Score")
        plt.close(fig)

        log_print("[Selesai] Seluruh analisis Tahap 3, 4, dan 5 serta Gambar PPI 03-06 telah disimpan.", log_fh)

if __name__ == "__main__":
    main()

