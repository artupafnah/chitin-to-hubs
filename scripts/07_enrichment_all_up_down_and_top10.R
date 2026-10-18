#!/usr/bin/env Rscript
# ==============================================================================
# Tahap 2 & Tahap 6:
# Analisis Pengayaan Fungsional (GO BP/MF/CC & KEGG) untuk Semua DEG, DEG Naik,
# DEG Turun, serta Top 10 Kandidat Hub PPI + Perbandingan Konsistensi
# Proyek: Analisis Transkriptomik GSE133053 (Kitin vs Kontrol Air)
# Dataset: NCBI GEO GSE133053 (Arabidopsis thaliana, Kitin C6 1h vs H2O 1h)
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  base_dir <- dirname(script_dir)
} else {
  base_dir <- getwd()
}

results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
logs_dir    <- file.path(base_dir, "logs")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "07_enrichment_all_up_down_and_top10.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, split = TRUE)

cat("==============================================================================\n")
cat("TAHAP 2 & TAHAP 6: ENRICHMENT GO (BP/MF/CC) & KEGG (ALL DEG, UP, DOWN, TOP 10)\n")
cat("==============================================================================\n")
cat("Tanggal Eksekusi:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# 1. Baca tabel lengkap DESeq2 (19.979 gen teruji sebagai Universe/Background)
deg_full_path <- file.path(results_dir, "DESeq2_all_genes_complete_table.csv")
df_all <- read.csv(deg_full_path, stringsAsFactors = FALSE)
universe_genes <- unique(df_all$gene_id)

df_deg_all  <- df_all %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) >= 1.0)
df_deg_up   <- df_all %>% filter(!is.na(padj), padj < 0.05, log2FC >= 1.0)
df_deg_down <- df_all %>% filter(!is.na(padj), padj < 0.05, log2FC <= -1.0)

genes_all  <- unique(df_deg_all$gene_id)
genes_up   <- unique(df_deg_up$gene_id)
genes_down <- unique(df_deg_down$gene_id)

# Baca Top 10 kandidat hub dari Tahap 5
top10_path <- file.path(results_dir, "top_candidates.csv")
df_top10 <- read.csv(top10_path, stringsAsFactors = FALSE)
genes_top10 <- unique(df_top10$gene_id)

# Pemetaan TAIR ID ke Simbol Gen agar kolom geneID lebih informatif
tair_to_sym <- setNames(df_all$symbol, df_all$gene_id)
convert_gene_ids_to_symbols <- function(gene_str) {
  ids <- unlist(strsplit(as.character(gene_str), "/"))
  syms <- ifelse(ids %in% names(tair_to_sym) & tair_to_sym[ids] != "", tair_to_sym[ids], ids)
  paste(syms, collapse = "/")
}

cat(sprintf("[Input] Universe (semua gen diuji): %d gen\n", length(universe_genes)))
cat(sprintf("[Input] Semua DEG (|log2FC|>=1.0, padj<0.05): %d gen\n", length(genes_all)))
cat(sprintf("[Input] DEG Naik (Up): %d gen | DEG Turun (Down): %d gen\n", length(genes_up), length(genes_down)))
cat(sprintf("[Input] Top 10 Kandidat Hub PPI: %d gen (%s)\n\n",
            length(genes_top10), paste(df_top10$symbol, collapse = ", ")))

# Fungsi bantu untuk menjalankan GO (BP, MF, CC) dan KEGG pada satu himpunan gen
run_full_enrichment <- function(gene_vec, set_label, universe_vec) {
  cat(sprintf("--- Menjalankan Enrichment untuk: %s (n = %d gen) ---\n", set_label, length(gene_vec)))
  out_list <- list()

  for (ont in c("BP", "MF", "CC")) {
    ego <- tryCatch(
      enrichGO(
        gene          = gene_vec,
        universe      = universe_vec,
        OrgDb         = org.At.tair.db,
        keyType       = "TAIR",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 1.0,
        qvalueCutoff  = 1.0,
        readable      = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      df_e <- as.data.frame(ego)
      df_e$Database <- paste0("GO:", ont)
      df_e$GeneSet  <- set_label
      df_e$geneSymbol <- vapply(df_e$geneID, convert_gene_ids_to_symbols, character(1))
      df_e$Significance <- ifelse(df_e$p.adjust < 0.05, "FDR < 0.05 (Signifikan)",
                           ifelse(df_e$pvalue < 0.05, "Nominal p < 0.05 (FDR >= 0.05)", "Tidak Signifikan"))
      n_fdr <- sum(df_e$p.adjust < 0.05)
      n_nom <- sum(df_e$pvalue < 0.05)
      cat(sprintf("   [GO:%s] Signifikan FDR < 0.05: %d term | Nominal p < 0.05: %d term\n", ont, n_fdr, n_nom))
      out_list[[paste0("GO_", ont)]] <- df_e
    } else {
      cat(sprintf("   [GO:%s] 0 term ditemukan\n", ont))
    }
  }

  # KEGG Enrichment (organism = 'ath')
  ekegg <- tryCatch(
    enrichKEGG(
      gene          = gene_vec,
      universe      = universe_vec,
      organism      = "ath",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1.0,
      qvalueCutoff  = 1.0
    ),
    error = function(e) NULL
  )
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    df_k <- as.data.frame(ekegg)
    df_k$Description <- sub(" - Arabidopsis thaliana \\(thale cress\\)$", "", df_k$Description)
    df_k$Database <- "KEGG"
    df_k$GeneSet  <- set_label
    df_k$geneSymbol <- vapply(df_k$geneID, convert_gene_ids_to_symbols, character(1))
    df_k$Significance <- ifelse(df_k$p.adjust < 0.05, "FDR < 0.05 (Signifikan)",
                         ifelse(df_k$pvalue < 0.05, "Nominal p < 0.05 (FDR >= 0.05)", "Tidak Signifikan"))
    # Samakan kolom bila ada tambahan category/subcategory di KEGG
    common_cols <- c("Database", "GeneSet", "ID", "Description", "GeneRatio", "BgRatio",
                     "pvalue", "p.adjust", "qvalue", "Count", "Significance", "geneID", "geneSymbol")
    for (cc in common_cols) if (!cc %in% names(df_k)) df_k[[cc]] <- NA
    df_k <- df_k[, common_cols]
    n_fdr_k <- sum(df_k$p.adjust < 0.05)
    n_nom_k <- sum(df_k$pvalue < 0.05)
    cat(sprintf("   [KEGG]  Signifikan FDR < 0.05: %d jalur | Nominal p < 0.05: %d jalur\n", n_fdr_k, n_nom_k))
    out_list[["KEGG"]] <- df_k
  } else {
    cat("   [KEGG]  0 jalur ditemukan\n")
  }

  common_cols <- c("Database", "GeneSet", "ID", "Description", "GeneRatio", "BgRatio",
                   "pvalue", "p.adjust", "qvalue", "Count", "Significance", "geneID", "geneSymbol")
  norm_list <- lapply(out_list, function(x) {
    for (cc in common_cols) if (!cc %in% names(x)) x[[cc]] <- NA
    x[, common_cols]
  })
  df_comb <- do.call(rbind, norm_list)
  rownames(df_comb) <- NULL
  return(df_comb)
}

res_all   <- run_full_enrichment(genes_all,   "All_DEG_717",   universe_genes)
res_up    <- run_full_enrichment(genes_up,    "Up_DEG_616",    universe_genes)
res_down  <- run_full_enrichment(genes_down,  "Down_DEG_101",  universe_genes)
res_top10 <- run_full_enrichment(genes_top10, "Top10_Hubs_10", universe_genes)

# Simpan tabel enrichment (simpan semua yang pvalue < 0.05 dengan kolom Significance eksplisit,
# serta tabel khusus FDR < 0.05)
save_enrichment_tables <- function(df_res, file_prefix) {
  df_sig_or_nom <- df_res %>%
    filter(pvalue < 0.05) %>%
    arrange(Database, p.adjust, pvalue)
  write.csv(df_sig_or_nom, file.path(results_dir, paste0(file_prefix, ".csv")), row.names = FALSE)
}

save_enrichment_tables(res_all,   "enrichment_ALL_DEG_GO_KEGG")
save_enrichment_tables(res_up,    "enrichment_UP_DEG_GO_KEGG")
save_enrichment_tables(res_down,  "enrichment_DOWN_DEG_GO_KEGG")
save_enrichment_tables(res_top10, "enrichment_TOP10_hubs_GO_KEGG")

# Buat tabel ringkasan jumlah term signifikan per kategori
all_sets_df <- bind_rows(res_all, res_up, res_down, res_top10)
summary_table <- all_sets_df %>%
  group_by(GeneSet, Database) %>%
  summarise(
    Total_Tested_Terms = n(),
    Significant_FDR_005 = sum(p.adjust < 0.05, na.rm = TRUE),
    Nominal_P_005 = sum(pvalue < 0.05, na.rm = TRUE),
    Top_Term_Description = Description[which.min(pvalue)][1],
    Top_Term_ID = ID[which.min(pvalue)][1],
    Top_Term_FDR = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_table, file.path(results_dir, "enrichment_summary_counts_by_category.csv"), row.names = FALSE)

cat("\n[Ringkasan Jumlah Term Enrichment Signifikan (FDR < 0.05 & Nominal p < 0.05)]:\n")
print(as.data.frame(summary_table))

# Buat tabel perbandingan langsung: Top 10 Hubs vs All 717 DEGs (Tahap 6)
top10_sig <- res_top10 %>% filter(p.adjust < 0.05 | (Database == "KEGG" & pvalue < 0.05))
all_lookup <- res_all %>%
  select(Database, ID, AllDEG_GeneRatio = GeneRatio, AllDEG_Count = Count,
         AllDEG_pvalue = pvalue, AllDEG_padj = p.adjust, AllDEG_Significance = Significance)

comp_top10_vs_all <- top10_sig %>%
  select(Database, ID, Description,
         Top10_GeneRatio = GeneRatio, Top10_Count = Count,
         Top10_pvalue = pvalue, Top10_padj = p.adjust,
         Top10_Significance = Significance, Top10_Genes = geneSymbol) %>%
  left_join(all_lookup, by = c("Database", "ID")) %>%
  mutate(
    Konsistensi_Status = case_when(
      !is.na(AllDEG_padj) & AllDEG_padj < 0.05 & Top10_padj < 0.05 ~ "Konsisten Signifikan di Keduanya (FDR < 0.05)",
      !is.na(AllDEG_padj) & AllDEG_padj >= 0.05 & Top10_padj < 0.05 ~ "Spesifik Signifikan pada Top 10 Hub (Terkonsentrasi di Hub)",
      TRUE ~ "Nominal / Diagnostik"
    )
  ) %>%
  arrange(Database, Top10_padj)

write.csv(comp_top10_vs_all, file.path(results_dir, "enrichment_comparison_Top10_vs_AllDEG.csv"), row.names = FALSE)

cat("\n[Tahap 6] Perbandingan Term Enrichment Top 10 Hub vs Semua DEG (717 gen):\n")
cat(sprintf("   - Total term signifikan (FDR < 0.05) pada Top 10 Hub: %d term\n", sum(res_top10$p.adjust < 0.05)))
cat(sprintf("   - Term yang juga signifikan (FDR < 0.05) pada Semua 717 DEG: %d term\n",
            sum(comp_top10_vs_all$Konsistensi_Status == "Konsisten Signifikan di Keduanya (FDR < 0.05)", na.rm = TRUE)))
cat(sprintf("   - Term yang spesifik mencapai FDR < 0.05 pada Top 10 Hub: %d term\n",
            sum(comp_top10_vs_all$Konsistensi_Status == "Spesifik Signifikan pada Top 10 Hub (Terkonsentrasi di Hub)", na.rm = TRUE)))

# Fungsi bantu parse GeneRatio (mis. "35/540" -> 0.0648)
parse_ratio <- function(r_chr) {
  vapply(strsplit(as.character(r_chr), "/"), function(x) {
    if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA_real_
  }, numeric(1))
}

save_ggplot_multi <- function(p, base_name, w = 15, h = 10) {
  ggsave(file.path(figures_dir, paste0(base_name, ".png")), plot = p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(figures_dir, paste0(base_name, ".pdf")), plot = p, width = w, height = h, bg = "white")
  ggsave(file.path(figures_dir, paste0(base_name, ".svg")), plot = p, width = w, height = h, bg = "white")
}

# ==============================================================================
# VISUALISASI 1: Fig_PPI_01_Enrichment_All_DEG_GO_KEGG
# Dotplot Top 10-12 Term per Kategori (GO:BP, GO:MF, GO:CC, KEGG) untuk Semua DEG (717 gen)
# ==============================================================================
prepare_top_terms <- function(df_sub, db_name, n_top = 12) {
  d <- df_sub %>% filter(Database == db_name)
  if (sum(d$p.adjust < 0.05) > 0) {
    d <- d %>% filter(p.adjust < 0.05) %>% arrange(p.adjust, pvalue) %>% head(n_top)
  } else {
    d <- d %>% filter(pvalue < 0.05) %>% arrange(pvalue) %>% head(n_top)
  }
  d$RatioNum <- parse_ratio(d$GeneRatio)
  d$NegLog10Padj <- -log10(pmax(d$p.adjust, 1e-300))
  d$ShortDesc <- ifelse(nchar(d$Description) > 48, paste0(substr(d$Description, 1, 45), "..."), d$Description)
  d$TermLabel <- paste0(d$ShortDesc, " (", d$ID, ")")
  return(d)
}

df_plot_all <- bind_rows(
  prepare_top_terms(res_all, "GO:BP", 12),
  prepare_top_terms(res_all, "GO:MF", 12),
  prepare_top_terms(res_all, "GO:CC", 10),
  prepare_top_terms(res_all, "KEGG",  10)
)
df_plot_all$DatabaseFacet <- factor(
  df_plot_all$Database,
  levels = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
  labels = c(
    "A. GO: Biological Process (FDR < 0.05)",
    "B. GO: Molecular Function (FDR < 0.05)",
    "C. GO: Cellular Component (Diagnostik Nominal p < 0.05; FDR >= 0.05)",
    "D. KEGG Pathway (2 Jalur Kanonik FDR < 0.05 + Top Nominal)"
  )
)

# Urutkan TermLabel berdasarkan RatioNum di tiap facet
df_plot_all <- df_plot_all %>%
  group_by(DatabaseFacet) %>%
  arrange(RatioNum, .by_group = TRUE) %>%
  mutate(TermLabel = factor(TermLabel, levels = unique(TermLabel))) %>%
  ungroup()

p_fig01 <- ggplot(df_plot_all, aes(x = RatioNum, y = TermLabel)) +
  geom_segment(aes(x = 0, xend = RatioNum, y = TermLabel, yend = TermLabel), color = "grey75", linewidth = 0.6) +
  geom_point(aes(size = Count, fill = NegLog10Padj, shape = Significance), color = "black", stroke = 0.6) +
  facet_wrap(~ DatabaseFacet, scales = "free", ncol = 2) +
  scale_fill_viridis_c(option = "plasma", name = expression(-log[10](p[adj]))) +
  scale_shape_manual(
    values = c("FDR < 0.05 (Signifikan)" = 21, "Nominal p < 0.05 (FDR >= 0.05)" = 24),
    name = "Status Signifikansi"
  ) +
  scale_size_continuous(range = c(3.2, 8.5), name = "Jumlah DEG") +
  labs(
    title = "Pengayaan Fungsional (GO:BP, GO:MF, GO:CC, & KEGG) Seluruh 717 DEG (Perlakuan Kitin C6 1h vs Kontrol H2O 1h)",
    subtitle = "Organisme: Arabidopsis thaliana (TAIR10) | Latar Belakang: 19.979 Gen Teruji | Koreksi Benjamini-Hochberg (FDR < 0.05)",
    x = "GeneRatio (Proporsi DEG terhadap Total Gen Teranotasi dalam Kategori)",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey25"),
    strip.background = element_rect(fill = "#1F4E78", color = "#1F4E78"),
    strip.text = element_text(color = "white", face = "bold", size = 10),
    axis.text.y = element_text(size = 8.8, color = "black"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.8))

save_ggplot_multi(p_fig01, "Fig_PPI_01_Enrichment_All_DEG_GO_KEGG", w = 16.5, h = 11)

# ==============================================================================
# VISUALISASI 2: Fig_PPI_02_Enrichment_Up_vs_Down_Comparison
# Perbandingan Enrichment Terpisah: DEG Naik (616 gen) vs DEG Turun (101 gen)
# ==============================================================================
df_up_top <- bind_rows(
  prepare_top_terms(res_up, "GO:BP", 10),
  prepare_top_terms(res_up, "GO:MF", 8),
  prepare_top_terms(res_up, "KEGG",  4)
) %>% mutate(GroupPanel = "A. DEG Naik pada Kitin C6 (616 Gen; GO:BP, GO:MF, & KEGG Signifikan FDR < 0.05)")

df_down_top <- bind_rows(
  prepare_top_terms(res_down, "GO:BP", 10),
  prepare_top_terms(res_down, "GO:MF", 6),
  prepare_top_terms(res_down, "GO:CC", 4)
) %>% mutate(GroupPanel = "B. DEG Turun pada Kitin C6 (101 Gen; Profil Diagnostik Nominal p < 0.05, FDR >= 0.05)")

df_up_down_plot <- bind_rows(df_up_top, df_down_top) %>%
  mutate(
    NegLog10P = -log10(pmax(pvalue, 1e-300)),
    TermWithDB = paste0("[", Database, "] ", ShortDesc, " (", ID, ")")
  ) %>%
  group_by(GroupPanel) %>%
  arrange(NegLog10P, .by_group = TRUE) %>%
  mutate(TermWithDB = factor(TermWithDB, levels = unique(TermWithDB))) %>%
  ungroup()

p_fig02 <- ggplot(df_up_down_plot, aes(x = NegLog10P, y = TermWithDB, fill = Database)) +
  geom_col(aes(alpha = Significance), color = "black", width = 0.72, linewidth = 0.35) +
  geom_text(aes(label = paste0("n=", Count, " (FDR=", formatC(p.adjust, format = "e", digits = 1), ")")),
            hjust = -0.05, size = 2.9, fontface = "bold") +
  facet_wrap(~ GroupPanel, scales = "free", ncol = 2) +
  scale_fill_manual(
    values = c("GO:BP" = "#D55E00", "GO:MF" = "#0072B2", "GO:CC" = "#009E73", "KEGG" = "#CC79A7"),
    name = "Kategori Database"
  ) +
  scale_alpha_manual(
    values = c("FDR < 0.05 (Signifikan)" = 0.95, "Nominal p < 0.05 (FDR >= 0.05)" = 0.55),
    name = "Status Signifikansi"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
  labs(
    title = "Perbandingan Pengayaan Fungsional Terpisah: DEG Naik (Up) vs DEG Turun (Down) pada Perlakuan Kitin C6",
    subtitle = "Arah Kontras: Perlakuan Kitin C6 1 jam vs Kontrol H2O 1 jam | Latar: 19.979 Gen Teruji",
    x = expression(-log[10](p-value~mentah)),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 10, color = "grey25"),
    strip.background = element_rect(fill = "#2F5597", color = "#2F5597"),
    strip.text = element_text(color = "white", face = "bold", size = 9.5),
    axis.text.y = element_text(size = 8.5, color = "black"),
    legend.position = "bottom"
  )

save_ggplot_multi(p_fig02, "Fig_PPI_02_Enrichment_Up_vs_Down_Comparison", w = 16.5, h = 9.5)

# ==============================================================================
# VISUALISASI 3: Fig_PPI_07_Top10_vs_AllDEG_Enrichment_Comparison
# Perbandingan Konsistensi Enrichment Top 10 Hub Prioritas vs Seluruh 717 DEG (Tahap 6)
# ==============================================================================
top_comp_terms <- comp_top10_vs_all %>%
  filter(Top10_padj < 0.05 | (Database == "KEGG" & Top10_pvalue < 0.05)) %>%
  group_by(Database) %>%
  slice_head(n = 10) %>%
  ungroup()

df_long_comp <- bind_rows(
  top_comp_terms %>%
    transmute(
      Database, ID, Description, Konsistensi_Status,
      Kelompok = "Top 10 Kandidat Hub PPI (n = 10 Gen)",
      Ratio = parse_ratio(Top10_GeneRatio),
      Count = Top10_Count,
      padj  = Top10_padj,
      Genes = Top10_Genes
    ),
  top_comp_terms %>%
    transmute(
      Database, ID, Description, Konsistensi_Status,
      Kelompok = "Seluruh DEG (n = 717 Gen)",
      Ratio = parse_ratio(AllDEG_GeneRatio),
      Count = AllDEG_Count,
      padj  = AllDEG_padj,
      Genes = ""
    )
) %>%
  mutate(
    NegLog10Padj = -log10(pmax(padj, 1e-300)),
    ShortDesc = ifelse(nchar(Description) > 46, paste0(substr(Description, 1, 43), "..."), Description),
    TermLabel = paste0("[", Database, "] ", ShortDesc, " (", ID, ")")
  )

term_order <- df_long_comp %>%
  filter(Kelompok == "Top 10 Kandidat Hub PPI (n = 10 Gen)") %>%
  arrange(NegLog10Padj) %>%
  pull(TermLabel) %>%
  unique()
df_long_comp$TermLabel <- factor(df_long_comp$TermLabel, levels = term_order)

p_fig07 <- ggplot(df_long_comp, aes(x = NegLog10Padj, y = TermLabel, fill = Kelompok)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "black", linewidth = 0.35) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  geom_text(
    aes(label = ifelse(Kelompok == "Top 10 Kandidat Hub PPI (n = 10 Gen)", paste0("Hub: ", Genes), paste0("All DEG: n=", Count))),
    position = position_dodge(width = 0.75),
    hjust = -0.04, size = 2.8, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Top 10 Kandidat Hub PPI (n = 10 Gen)" = "#D55E00", "Seluruh DEG (n = 717 Gen)" = "#0072B2"),
    name = "Himpunan Gen"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.42))) +
  labs(
    title = "Perbandingan Pengayaan Fungsional (GO & KEGG): Top 10 Kandidat Hub PPI vs Seluruh 717 DEG",
    subtitle = "Garis putus-putus merah menunjukkan batas signifikansi FDR = 0.05 (-log10(padj) = 1.30) | Kontras: Kitin C6 1h vs Kontrol H2O 1h",
    x = expression(-log[10](p[adj]~Benjamini-Hochberg)),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5),
    plot.subtitle = element_text(size = 9.8, color = "grey25"),
    axis.text.y = element_text(size = 8.8, color = "black"),
    legend.position = "bottom"
  )

save_ggplot_multi(p_fig07, "Fig_PPI_07_Top10_vs_AllDEG_Enrichment_Comparison", w = 15.5, h = 9.2)

cat("\n[Selesai] Seluruh tabel enrichment (Tahap 2 & 6) dan Gambar PPI 01, 02, 07 telah disimpan.\n")
sink()
close(log_con)

