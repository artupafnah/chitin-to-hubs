#!/usr/bin/env Rscript
# ==============================================================================
# TAHAP 4 (Visualisasi QC), TAHAP 5 (Filtering & Normalisasi),
# TAHAP 6 (Eksplorasi Sampel), TAHAP 7 (Uji DEG DESeq2 + apeglm),
# DAN TAHAP 8 (Visualisasi Publikasi 300 DPI PNG + PDF + SVG)
#
# Dataset   : GSE133053 / BioProject PRJNA549850 (Wan et al., 2019, Nature Plants)
# Organisme : Arabidopsis thaliana (Col-0, TAIR10)
# Kontras   : Perlakuan Elisitor Kitin 1 jam (C6) vs Kontrol Air 1 jam (H2O)
#             Arah kontras SELALU Perlakuan (C6) vs Kontrol (H2O)
#             (log2FC positif = ekspresi NAIK pada perlakuan elisitor kitin)
# Seed      : 42
# ==============================================================================

set.seed(42)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(grid)
  library(dplyr)
  library(tidyr)
})

cat("[INFO] Memulai analisis Tahap 4 s.d. Tahap 8 (DESeq2 RNA-seq GSE133053)...\n")

# Palet warna ramah buta warna (Okabe-Ito)
col_up   <- "#D55E00" # Vermilion (Naik pada Perlakuan Kitin C6)
col_down <- "#0072B2" # Blue (Turun pada Perlakuan Kitin C6)
col_ns   <- "#999999" # Grey (Tidak signifikan)
col_c6   <- "#E69F00" # Orange (Perlakuan Kitin C6)
col_h2o  <- "#56B4E9" # Sky Blue (Kontrol H2O)

save_pub_plot <- function(plot_obj, base_name, width = 9, height = 6) {
  png_file <- file.path("figures", paste0(base_name, ".png"))
  pdf_file <- file.path("figures", paste0(base_name, ".pdf"))
  svg_file <- file.path("figures", paste0(base_name, ".svg"))
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_file, plot = plot_obj, width = width, height = height, bg = "white")
  tryCatch({
    ggsave(svg_file, plot = plot_obj, width = width, height = height, bg = "white")
  }, error = function(e) {
    # Fallback jika svglite tidak tersedia: gunakan perangkat svg bawaan R
    svg(svg_file, width = width, height = height)
    print(plot_obj)
    dev.off()
  })
  cat(sprintf("[GAMBAR] Tersimpan: %s (.png 300 dpi, .pdf, .svg)\n", base_name))
}

# ==============================================================================
# 1. MUAT METADATA, TABEL QC, DAN MATRIKS COUNT (TAHAP 4)
# ==============================================================================
qc_df   <- read.csv("qc/qc_library_and_mapping_summary.csv")
meta_8  <- read.csv("data/metadata.csv")
meta_7  <- read.csv("data/metadata_valid7samples.csv")

raw_8_df <- read.csv("results/raw_counts_matrix_all8samples.csv", row.names = 1, check.names = FALSE)
raw_7_df <- read.csv("results/raw_counts_matrix_valid7samples.csv", row.names = 1, check.names = FALSE)

mat_8 <- as.matrix(raw_8_df)
mat_7 <- as.matrix(raw_7_df)

# --- GAMBAR 1: Visualisasi QC Kedalaman Sekuensing & Persentase Mapping ---
qc_df$sample_id <- factor(qc_df$sample_id, levels = qc_df$sample_id)
qc_df$Kondisi_Label <- ifelse(qc_df$kondisi == "Perlakuan_Kitin_C6",
                              "Perlakuan: Kitin C6 (1 jam)", "Kontrol: Air H2O (1 jam)")

p1_a <- ggplot(qc_df, aes(x = sample_id, y = assigned_reads_sum_isoforms / 1e6, fill = Kondisi_Label)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f M", assigned_reads_sum_isoforms / 1e6)),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("Kontrol: Air H2O (1 jam)" = col_h2o,
                               "Perlakuan: Kitin C6 (1 jam)" = col_c6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "A. Jumlah Read Ter-assign Unik per Sampel (TAIR10)",
    subtitle = "Kontras: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam)",
    x = NULL, y = "Read Ter-assign (Juta / Millions)", fill = "Grup Eksperimen"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
        legend.position = "top",
        plot.title = element_text(face = "bold"))

p1_b <- ggplot(qc_df, aes(x = sample_id, y = pct_assigned_sum, fill = qc_flag_threshold_70pct)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  annotate("text", x = 4.5, y = 74, label = "Ambang Peringatan QC (70%)",
           color = "#D55E00", fontface = "bold", size = 3.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct_assigned_sum)),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("OK (>=70%)" = "#009E73", "WARNING (<70%)" = "#D55E00"),
                    name = "Status Mapping QC") +
  scale_y_continuous(limits = c(0, 108)) +
  labs(
    title = "B. Persentase Read Ter-assign terhadap Total Read SRA",
    subtitle = "C6_1h_Rep1 (0,64%; 4.816 read) gagal sekuensing (dieksklusi); C6_1h_Rep3 (12,29%; 2,38 juta read) dipertahankan",
    x = "ID Sampel", y = "Persentase Mapping / Assignment (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
        legend.position = "top",
        plot.title = element_text(face = "bold"))

fig1 <- p1_a / p1_b + plot_annotation(
  title = "Gambar 1. Evaluasi Quality Control (QC) Kuantifikasi Read RNA-seq (GSE133053)",
  theme = theme(plot.title = element_text(face = "bold", size = 13))
)
save_pub_plot(fig1, "Fig01_QC_Library_and_Mapping_Rate", width = 10, height = 8)

# ==============================================================================
# 2. TAHAP 5: FILTERING GEN RENDAH DAN NORMALISASI DESeq2
# ==============================================================================
rownames(meta_8) <- meta_8$sample_id
meta_8$kondisi <- factor(meta_8$kondisi, levels = c("Kontrol_H2O", "Perlakuan_Kitin_C6"))

rownames(meta_7) <- meta_7$sample_id
meta_7$kondisi <- factor(meta_7$kondisi, levels = c("Kontrol_H2O", "Perlakuan_Kitin_C6"))
meta_7$lane_id <- factor(meta_7$lane_id)

# Buat objek DESeqDataSet untuk 8 sampel (diagnostik QC) dan 7 sampel valid (analisis utama)
dds_8_raw <- DESeqDataSetFromMatrix(countData = mat_8, colData = meta_8, design = ~ kondisi)
dds_7_raw <- DESeqDataSetFromMatrix(countData = mat_7, colData = meta_7, design = ~ kondisi)

# Kriteria filtering: count >= 10 pada minimal n sampel (n = 3, ukuran grup terkecil pada C6)
min_samples <- 3
keep_genes_7 <- rowSums(counts(dds_7_raw) >= 10) >= min_samples
dds_7 <- dds_7_raw[keep_genes_7, ]

keep_genes_8 <- rowSums(counts(dds_8_raw) >= 10) >= min_samples
dds_8 <- dds_8_raw[keep_genes_8, ]

cat(sprintf("[TAHAP 5] Total lokus gen awal TAIR10 : %d gen\n", nrow(dds_7_raw)))
cat(sprintf("[TAHAP 5] Gen lolos filtering (>=10 di >=3 sampel) : %d gen (Dibuang: %d gen)\n",
            nrow(dds_7), nrow(dds_7_raw) - nrow(dds_7)))

# Estimasi size factor (median-of-ratios) & VST
dds_8 <- estimateSizeFactors(dds_8)
dds_7 <- estimateSizeFactors(dds_7)

vsd_8 <- vst(dds_8, blind = TRUE)
vsd_7_blind <- vst(dds_7, blind = TRUE)
vsd_7 <- vst(dds_7, blind = FALSE)

norm_counts_7 <- counts(dds_7, normalized = TRUE)
vst_mat_7     <- assay(vsd_7)

write.csv(norm_counts_7, "results/normalized_counts_DESeq2.csv")
write.csv(vst_mat_7, "results/vst_transformed_matrix.csv")

sf_summary <- data.frame(
  sample_id = colnames(dds_7),
  kondisi = as.character(meta_7$kondisi),
  raw_library_size_all_genes = colSums(mat_7),
  raw_library_size_filtered_genes = colSums(counts(dds_7)),
  DESeq2_size_factor = round(sizeFactors(dds_7), 5),
  normalized_library_size = round(colSums(norm_counts_7), 1)
)
write.csv(sf_summary, "results/size_factors_and_filtering_summary.csv", row.names = FALSE)
print(sf_summary)

# --- GAMBAR 2: Boxplot Distribusi Ekspresi SEBELUM vs SESUDAH Normalisasi (Satu Panel) ---
df_before <- as.data.frame(log2(counts(dds_7) + 1)) %>%
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "expr") %>%
  mutate(Stage = "Sebelum Normalisasi: log2(Raw Count + 1)")

df_after_norm <- as.data.frame(log2(norm_counts_7 + 1)) %>%
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "expr") %>%
  mutate(Stage = "Sesudah Normalisasi: log2(DESeq2 Normalized Count + 1)")

df_after_vst <- as.data.frame(vst_mat_7) %>%
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "expr") %>%
  mutate(Stage = "Sesudah Transformasi: Variance Stabilizing Transformation (VST)")

df_box <- bind_rows(df_before, df_after_norm, df_after_vst) %>%
  left_join(meta_7[, c("sample_id", "kondisi")], by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = meta_7$sample_id),
    Stage = factor(Stage, levels = c(
      "Sebelum Normalisasi: log2(Raw Count + 1)",
      "Sesudah Normalisasi: log2(DESeq2 Normalized Count + 1)",
      "Sesudah Transformasi: Variance Stabilizing Transformation (VST)"
    )),
    Kondisi_Label = ifelse(kondisi == "Perlakuan_Kitin_C6",
                           "Perlakuan: Kitin C6 (1 jam)", "Kontrol: Air H2O (1 jam)")
  )

fig2 <- ggplot(df_box, aes(x = sample_id, y = expr, fill = Kondisi_Label)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25, width = 0.65) +
  facet_wrap(~ Stage, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("Kontrol: Air H2O (1 jam)" = col_h2o,
                               "Perlakuan: Kitin C6 (1 jam)" = col_c6)) +
  labs(
    title = "Gambar 2. Distribusi Ekspresi Gen Sebelum dan Sesudah Normalisasi DESeq2 (19.979 Gen Terfilter)",
    subtitle = "Kontras: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam) | Normalisasi Median-of-Ratios menyelaraskan median C6_1h_Rep3",
    x = "Sampel Biologis", y = "Tingkat Ekspresi (Skala Log2 / VST)", fill = "Kondisi Eksperimen"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "#F0F4F8"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12.5)
  )
save_pub_plot(fig2, "Fig02_Boxplot_Before_vs_After_Normalization", width = 13, height = 5.8)

# ==============================================================================
# 3. TAHAP 6: EKSPLORASI SAMPEL (PCA TOP 500 GEN & HEATMAP KORELASI / JARAK)
# ==============================================================================
compute_pca_df <- function(vsd_obj, meta_df, ntop = 500) {
  mat <- assay(vsd_obj)
  rv  <- rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(mat[select, ]), center = TRUE, scale. = FALSE)
  percentVar <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    kondisi = meta_df$kondisi,
    lane_id = meta_df$lane_id,
    Kondisi_Label = ifelse(meta_df$kondisi == "Perlakuan_Kitin_C6",
                           "Perlakuan: Kitin C6 (1 jam)", "Kontrol: Air H2O (1 jam)")
  )
  list(df = df, pct = percentVar)
}

pca_8 <- compute_pca_df(vsd_8, meta_8, ntop = 500)
pca_7 <- compute_pca_df(vsd_7_blind, meta_7, ntop = 500)

write.csv(pca_7$df, "results/PCA_coordinates_valid7samples.csv", row.names = FALSE)

p3_a <- ggplot(pca_8$df, aes(x = PC1, y = PC2, color = Kondisi_Label, shape = lane_id)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3.2, fontface = "bold", max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("Kontrol: Air H2O (1 jam)" = "#0072B2",
                                "Perlakuan: Kitin C6 (1 jam)" = "#D55E00")) +
  labs(
    title = "A. PCA Seluruh 8 Sampel Awal (Sebelum Eksklusi Sampel Gagal)",
    subtitle = "Outlier ekstrem C6_1h_Rep1 (hanya 4.816 read) mendominasi PC1",
    x = sprintf("PC1: %.1f%% Varians", pca_8$pct[1]),
    y = sprintf("PC2: %.1f%% Varians", pca_8$pct[2]),
    color = "Kondisi", shape = "Sequencing Lane"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p3_b <- ggplot(pca_7$df, aes(x = PC1, y = PC2, color = Kondisi_Label, shape = lane_id)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
  geom_point(size = 4.5, alpha = 0.95) +
  geom_text_repel(aes(label = sample_id), size = 3.4, fontface = "bold", box.padding = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c("Kontrol: Air H2O (1 jam)" = "#0072B2",
                                "Perlakuan: Kitin C6 (1 jam)" = "#D55E00")) +
  labs(
    title = "B. PCA 7 Sampel Valid Pasca-QC (Top 500 Gen Paling Variabel)",
    subtitle = "PC1 memisahkan secara tegas Perlakuan Kitin C6 (kanan) vs Kontrol H2O (kiri)",
    x = sprintf("PC1: %.1f%% Varians", pca_7$pct[1]),
    y = sprintf("PC2: %.1f%% Varians", pca_7$pct[2]),
    color = "Kondisi", shape = "Sequencing Lane"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

fig3 <- (p3_a | p3_b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig3 <- fig3 + plot_annotation(
  title = "Gambar 3. Analisis Komponen Utama (PCA) Data Transformasi VST (Top 500 Gen Variabel)",
  subtitle = "Kontras: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam) | Dataset RNA-seq GSE133053",
  theme = theme(plot.title = element_text(face = "bold", size = 13))
)
save_pub_plot(fig3, "Fig03_PCA_Plot_VST_Top500", width = 12.5, height = 6)

# --- GAMBAR 4: Heatmap Korelasi Pearson & Jarak Euclidean Antar Sampel ---
cor_mat_7  <- cor(assay(vsd_7_blind), method = "pearson")
dist_mat_7 <- as.matrix(dist(t(assay(vsd_7_blind)), method = "euclidean"))

write.csv(round(cor_mat_7, 5), "results/sample_correlation_matrix_VST.csv")
write.csv(round(dist_mat_7, 3), "results/sample_euclidean_distance_VST.csv")

ann_col_sample <- data.frame(
  Kondisi = ifelse(meta_7$kondisi == "Perlakuan_Kitin_C6", "Perlakuan Kitin (C6)", "Kontrol Air (H2O)"),
  Lane = as.character(meta_7$lane_id),
  row.names = rownames(meta_7)
)
ann_colors_sample <- list(
  Kondisi = c("Kontrol Air (H2O)" = col_h2o, "Perlakuan Kitin (C6)" = col_c6),
  Lane = c("Lane_1" = "#CC79A7", "Lane_2" = "#009E73", "Lane_3" = "#F0E442", "Lane_4" = "#0072B2")
)

save_pheatmap_all <- function(expr_call, base_name, width = 8, height = 6.5) {
  png(file.path("figures", paste0(base_name, ".png")), width = width * 300, height = height * 300, res = 300)
  eval(expr_call)
  dev.off()
  pdf(file.path("figures", paste0(base_name, ".pdf")), width = width, height = height)
  eval(expr_call)
  dev.off()
  svg(file.path("figures", paste0(base_name, ".svg")), width = width, height = height)
  eval(expr_call)
  dev.off()
  cat(sprintf("[GAMBAR] Tersimpan: %s (.png 300 dpi, .pdf, .svg)\n", base_name))
}

cor_call <- quote(pheatmap(
  cor_mat_7,
  annotation_col = ann_col_sample,
  annotation_row = ann_col_sample[, "Kondisi", drop = FALSE],
  annotation_colors = ann_colors_sample,
  display_numbers = TRUE,
  number_format = "%.3f",
  number_color = "black",
  fontsize_number = 9,
  color = viridis(100, option = "viridis"),
  clustering_distance_rows = as.dist(1 - cor_mat_7),
  clustering_distance_cols = as.dist(1 - cor_mat_7),
  main = "Gambar 4. Heatmap Korelasi Pearson Antar Sampel (VST, Klaster Hierarkis)\nKontras: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam)"
))
save_pheatmap_all(cor_call, "Fig04_Sample_Correlation_Distance_Heatmap", width = 8.5, height = 6.8)

# ==============================================================================
# 4. TAHAP 7: UJI DIFERENSIAL EKSPRESI GEN (DESeq2 + Wald Test + BH + apeglm)
# ==============================================================================
cat("\n[TAHAP 7] Menjalankan DESeq2 Wald test dan lfcShrink (apeglm)...\n")
set.seed(42)
dds_7 <- DESeq(dds_7)

# Pastikan nama koefisien kontras adalah kondisi_Perlakuan_Kitin_C6_vs_Kontrol_H2O
cat("[INFO] Nama koefisien model DESeq2:", resultsNames(dds_7), "\n")
res_unshrunk <- results(dds_7, contrast = c("kondisi", "Perlakuan_Kitin_C6", "Kontrol_H2O"), alpha = 0.05)
res_shrunk   <- lfcShrink(dds_7, coef = "kondisi_Perlakuan_Kitin_C6_vs_Kontrol_H2O",
                          res = res_unshrunk, type = "apeglm")

# Buat tabel hasil lengkap dengan anotasi simbol & deskripsi gen dari org.At.tair.db
gene_ids <- rownames(res_shrunk)
symbols <- mapIds(org.At.tair.db, keys = gene_ids, column = "SYMBOL", keytype = "TAIR", multiVals = "first")
genenames <- mapIds(org.At.tair.db, keys = gene_ids, column = "GENENAME", keytype = "TAIR", multiVals = "first")

res_df <- data.frame(
  gene_id = gene_ids,
  symbol = ifelse(is.na(symbols) | symbols == "", gene_ids, as.character(symbols)),
  description = ifelse(is.na(genenames) | genenames == "", "Uncharacterized Arabidopsis protein (TAIR10)", as.character(genenames)),
  baseMean = round(res_shrunk$baseMean, 3),
  log2FC = round(res_shrunk$log2FoldChange, 4),          # Shrunk log2FC (apeglm) - UTAMA
  lfcSE = round(res_shrunk$lfcSE, 4),                    # Posterior SD dari apeglm
  log2FC_unshrunk = round(res_unshrunk$log2FoldChange, 4),
  lfcSE_unshrunk = round(res_unshrunk$lfcSE, 4),
  stat = round(res_unshrunk$stat, 4),                    # Wald statistic Z
  pvalue = res_unshrunk$pvalue,
  padj = res_unshrunk$padj,
  contrast_direction = "Perlakuan_Kitin_C6_1h_vs_Kontrol_H2O_1h (log2FC > 0 = Naik pada Kitin)",
  stringsAsFactors = FALSE
)

# Kategori regulasi pada ambang utama: padj < 0.05 & |log2FC| >= 1.0 (shrunk apeglm)
res_df$regulation <- "Tidak Signifikan"
res_df$regulation[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FC >= 1.0]  <- "Naik pada Perlakuan Kitin (Up)"
res_df$regulation[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FC <= -1.0] <- "Turun pada Perlakuan Kitin (Down)"

# Urutkan dari padj terkecil lalu |log2FC| terbesar
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FC)), ]

# Simpan tabel lengkap dan daftar DEG
deg_all  <- subset(res_df, regulation != "Tidak Signifikan")
deg_up   <- subset(res_df, regulation == "Naik pada Perlakuan Kitin (Up)")
deg_down <- subset(res_df, regulation == "Turun pada Perlakuan Kitin (Down)")
deg_up   <- deg_up[order(-deg_up$log2FC), ]
deg_down <- deg_down[order(deg_down$log2FC), ]

write.csv(res_df,   "results/DESeq2_all_genes_complete_table.csv", row.names = FALSE)
write.csv(deg_all,  "results/DEG_significant_all_logFC1.csv", row.names = FALSE)
write.csv(deg_up,   "results/DEG_upregulated_in_Chitin_C6.csv", row.names = FALSE)
write.csv(deg_down, "results/DEG_downregulated_in_Chitin_C6.csv", row.names = FALSE)

cat(sprintf("\n================ RINGKASAN HASIL DEG UTAMA (padj < 0.05 & |log2FC| >= 1.0) ================\n"))
cat(sprintf("Total Gen Diuji                : %d gen\n", nrow(res_df)))
cat(sprintf("Total DEG Signifikan (apeglm)  : %d gen\n", nrow(deg_all)))
cat(sprintf("  - NAIK pada Kitin (log2FC>=1): %d gen\n", nrow(deg_up)))
cat(sprintf("  - TURUN pada Kitin (log2FC<=-1): %d gen\n", nrow(deg_down)))

# Analisis Sensitivitas pada |log2FC| >= 0.58, 1.0, dan 1.5 (Shrunk apeglm vs Unshrunk MLE)
# Serta perbandingan model 8 sampel vs 7 sampel valid vs 6 sampel (tanpa C6_1h_Rep3)
dds_8_run <- DESeq(dds_8, quiet = TRUE)
res_8_un  <- results(dds_8_run, contrast = c("kondisi", "Perlakuan_Kitin_C6", "Kontrol_H2O"), alpha = 0.05)
res_8_sh  <- lfcShrink(dds_8_run, coef = "kondisi_Perlakuan_Kitin_C6_vs_Kontrol_H2O", res = res_8_un, type = "apeglm", quiet = TRUE)

idx_6 <- colnames(mat_7) != "C6_1h_Rep3"
meta_6 <- meta_7[idx_6, ]
dds_6 <- DESeqDataSetFromMatrix(mat_7[, idx_6], colData = meta_6, design = ~ kondisi)
dds_6 <- dds_6[rowSums(counts(dds_6) >= 10) >= 2, ]
dds_6 <- DESeq(dds_6, quiet = TRUE)
res_6_un <- results(dds_6, contrast = c("kondisi", "Perlakuan_Kitin_C6", "Kontrol_H2O"), alpha = 0.05)
res_6_sh <- lfcShrink(dds_6, coef = "kondisi_Perlakuan_Kitin_C6_vs_Kontrol_H2O", res = res_6_un, type = "apeglm", quiet = TRUE)

count_deg <- function(r_obj, lfc_thr, model_label, method_label) {
  up <- sum(!is.na(r_obj$padj) & r_obj$padj < 0.05 & r_obj$log2FoldChange >= lfc_thr)
  dn <- sum(!is.na(r_obj$padj) & r_obj$padj < 0.05 & r_obj$log2FoldChange <= -lfc_thr)
  data.frame(
    Model_Set = model_label,
    Estimator_log2FC = method_label,
    Threshold_abs_log2FC = lfc_thr,
    Fold_Change_Linear = round(2^lfc_thr, 2),
    padj_cutoff = 0.05,
    DEG_Up_in_Chitin = up,
    DEG_Down_in_Chitin = dn,
    DEG_Total = up + dn
  )
}

sens_list <- list()
for (thr in c(0.58, 1.0, 1.5)) {
  sens_list[[length(sens_list) + 1]] <- count_deg(res_shrunk,   thr, "7 Sampel Valid (Utama: 4 H2O vs 3 C6)", "Shrunk (apeglm)")
  sens_list[[length(sens_list) + 1]] <- count_deg(res_unshrunk, thr, "7 Sampel Valid (Utama: 4 H2O vs 3 C6)", "Unshrunk (MLE)")
  sens_list[[length(sens_list) + 1]] <- count_deg(res_8_sh,     thr, "8 Sampel Awal (Termasuk C6_1h_Rep1 Gagal)", "Shrunk (apeglm)")
  sens_list[[length(sens_list) + 1]] <- count_deg(res_6_sh,     thr, "6 Sampel Ketat (Tanpa C6_1h_Rep1 & Rep3)", "Shrunk (apeglm)")
}
sens_df <- bind_rows(sens_list)
write.csv(sens_df, "results/DEG_sensitivity_summary_table.csv", row.names = FALSE)
cat("\n=== Tabel Sensitivitas Ambang |log2FC| dan Set Sampel ===\n")
print(sens_df)

# ==============================================================================
# 5. TAHAP 8: VISUALISASI HASIL DEG
# ==============================================================================

# --- GAMBAR 5: Volcano Plot (Tandai Naik/Turun, Label ~15 Gen Teratas, Garis Ambang) ---
volc_df <- res_df %>% filter(!is.na(padj) & !is.na(log2FC))
volc_df$neg_log10_padj <- -log10(pmax(volc_df$padj, 1e-300))

# Pilih 15 gen teratas yang paling signifikan & biologis menonjol untuk diberi label
top_up_labels   <- volc_df %>% filter(regulation == "Naik pada Perlakuan Kitin (Up)") %>%
  arrange(padj, desc(log2FC)) %>% head(10)
top_down_labels <- volc_df %>% filter(regulation == "Turun pada Perlakuan Kitin (Down)") %>%
  arrange(padj, log2FC) %>% head(5)
label_df <- bind_rows(top_up_labels, top_down_labels) %>%
  mutate(display_label = ifelse(symbol == gene_id, gene_id, paste0(symbol, " (", gene_id, ")")))

fig5 <- ggplot(volc_df, aes(x = log2FC, y = neg_log10_padj, color = regulation)) +
  geom_point(alpha = 0.65, size = 1.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30", linewidth = 0.6) +
  geom_label_repel(
    data = label_df,
    aes(label = display_label),
    size = 2.9, fontface = "bold", box.padding = 0.35, max.overlaps = 25,
    fill = "white", alpha = 0.9, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Naik pada Perlakuan Kitin (Up)" = col_up,
      "Turun pada Perlakuan Kitin (Down)" = col_down,
      "Tidak Signifikan" = col_ns
    ),
    labels = c(
      "Naik pada Perlakuan Kitin (Up)" = sprintf("Naik pada Kitin (n = %d)", nrow(deg_up)),
      "Turun pada Perlakuan Kitin (Down)" = sprintf("Turun pada Kitin (n = %d)", nrow(deg_down)),
      "Tidak Signifikan" = sprintf("Tidak Signifikan (n = %d)", nrow(volc_df) - nrow(deg_all))
    )
  ) +
  labs(
    title = "Gambar 5. Volcano Plot Ekspresi Gen Diferensial: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam)",
    subtitle = sprintf("DESeq2 Wald Test + apeglm Shrinkage | Kriteria: padj < 0,05 & |log2FC| >= 1,0 (Naik: %d gen, Turun: %d gen)",
                       nrow(deg_up), nrow(deg_down)),
    x = "log2 Fold Change (Perlakuan Kitin C6 vs Kontrol H2O; Positif = Naik pada Kitin)",
    y = "-log10(Adjusted P-value / FDR Benjamini-Hochberg)",
    color = "Status Regulasi DEG"
  ) +
  theme_bw(base_size = 11.5) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 12.5))
save_pub_plot(fig5, "Fig05_Volcano_Plot_C6_vs_H2O", width = 10, height = 7)

# --- GAMBAR 6: MA Plot (log2FC vs Rata-Rata Ekspresi, Sebelum vs Sesudah apeglm Shrinkage) ---
ma_df <- volc_df %>%
  mutate(log10_baseMean = log10(baseMean + 1))

p6_a <- ggplot(ma_df, aes(x = log10_baseMean, y = log2FC_unshrunk, color = regulation)) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"), color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Naik pada Perlakuan Kitin (Up)" = col_up,
                                "Turun pada Perlakuan Kitin (Down)" = col_down,
                                "Tidak Signifikan" = col_ns)) +
  coord_cartesian(ylim = c(-8, 10)) +
  labs(
    title = "A. Sebelum Shrinkage (Unshrunk MLE log2FC)",
    subtitle = "Gen ber-count rendah (kiri) memiliki varians log2FC ekstrem",
    x = "log10(Rata-rata Ekspresi Ternormalisasi / baseMean + 1)",
    y = "log2 Fold Change (Unshrunk MLE)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

p6_b <- ggplot(ma_df, aes(x = log10_baseMean, y = log2FC, color = regulation)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"), color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Naik pada Perlakuan Kitin (Up)" = col_up,
                                "Turun pada Perlakuan Kitin (Down)" = col_down,
                                "Tidak Signifikan" = col_ns)) +
  coord_cartesian(ylim = c(-8, 10)) +
  labs(
    title = "B. Sesudah Shrinkage Empirical Bayes (apeglm)",
    subtitle = "Noise pada baseMean rendah berhasil diredam sementara sinyal DEG sejati terjaga",
    x = "log10(Rata-rata Ekspresi Ternormalisasi / baseMean + 1)",
    y = "log2 Fold Change (Shrunk apeglm)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

fig6 <- (p6_a | p6_b) + plot_annotation(
  title = "Gambar 6. MA Plot (log2FC vs Mean of Normalized Counts): Perlakuan Kitin C6 (1 jam) vs Kontrol H2O (1 jam)",
  theme = theme(plot.title = element_text(face = "bold", size = 12.5))
)
save_pub_plot(fig6, "Fig06_MA_Plot_Before_and_After_apeglm", width = 12.5, height = 5.8)

# --- GAMBAR 7: Heatmap Top 50 DEG (Z-score per Gen, Anotasi Kolom Kondisi) ---
top50_deg <- deg_all %>% arrange(padj, desc(abs(log2FC))) %>% head(50)
mat_top50 <- vst_mat_7[top50_deg$gene_id, ]
# Z-score per baris (gen)
mat_top50_z <- t(scale(t(mat_top50)))

row_labels_50 <- ifelse(
  top50_deg$symbol == top50_deg$gene_id,
  top50_deg$gene_id,
  paste0(top50_deg$gene_id, " (", top50_deg$symbol, ")")
)
rownames(mat_top50_z) <- row_labels_50

ann_row_50 <- data.frame(
  Arah_DEG = ifelse(top50_deg$log2FC > 0, "Naik (Up in C6)", "Turun (Down in C6)"),
  row.names = row_labels_50
)
ann_colors_50 <- list(
  Kondisi = c("Kontrol Air (H2O)" = col_h2o, "Perlakuan Kitin (C6)" = col_c6),
  Lane = c("Lane_1" = "#CC79A7", "Lane_2" = "#009E73", "Lane_3" = "#F0E442", "Lane_4" = "#0072B2"),
  Arah_DEG = c("Naik (Up in C6)" = col_up, "Turun (Down in C6)" = col_down)
)

hm50_call <- quote(pheatmap(
  mat_top50_z,
  annotation_col = ann_col_sample,
  annotation_row = ann_row_50,
  annotation_colors = ann_colors_50,
  color = colorRampPalette(c("#0072B2", "white", "#D55E00"))(100),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7.5,
  fontsize_col = 10,
  border_color = "grey85",
  main = "Gambar 7. Heatmap Top 50 DEG Signifikan (Z-score Ekspresi VST per Gen)\nKontras: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam)"
))
save_pheatmap_all(hm50_call, "Fig07_Heatmap_Top50_DEG_Zscore", width = 9.5, height = 10.5)

# --- GAMBAR 8: Histogram Distribusi P-value & Adjusted P-value ---
pval_df <- volc_df %>%
  dplyr::select(pvalue, padj) %>%
  pivot_longer(cols = c(pvalue, padj), names_to = "Tipe", values_to = "val") %>%
  mutate(Tipe_Label = ifelse(Tipe == "pvalue",
                             "A. Raw P-value (Wald Test)",
                             "B. Adjusted P-value / FDR (Benjamini-Hochberg)"))

fig8 <- ggplot(pval_df, aes(x = val, fill = Tipe_Label)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.25, boundary = 0) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  facet_wrap(~ Tipe_Label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("A. Raw P-value (Wald Test)" = "#0072B2",
                               "B. Adjusted P-value / FDR (Benjamini-Hochberg)" = "#009E73")) +
  labs(
    title = "Gambar 8. Histogram Distribusi Raw P-value dan Adjusted P-value (19.979 Gen Teruji)",
    subtitle = "Puncak tajam mendekati 0 (anti-konservatif sehat) membuktikan adanya sinyal ekspresi diferensial sejati akibat perlakuan kitin",
    x = "Nilai Probabilitas (P-value / padj; Garis Putus-putus = Ambang 0,05)",
    y = "Frekuensi (Jumlah Gen)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#F0F4F8"),
        strip.text = element_text(face = "bold", size = 10.5),
        plot.title = element_text(face = "bold", size = 12.5))
save_pub_plot(fig8, "Fig08_Pvalue_Distribution_Histogram", width = 10.5, height = 5)

# --- GAMBAR 9: Bar Chart Jumlah DEG Naik vs Turun dan Analisis Sensitivitas Ambang |log2FC| ---
sens_main <- sens_df %>%
  filter(Model_Set == "7 Sampel Valid (Utama: 4 H2O vs 3 C6)" & Estimator_log2FC == "Shrunk (apeglm)") %>%
  dplyr::select(Threshold_abs_log2FC, DEG_Up_in_Chitin, DEG_Down_in_Chitin) %>%
  pivot_longer(cols = c(DEG_Up_in_Chitin, DEG_Down_in_Chitin), names_to = "Arah", values_to = "Jumlah") %>%
  mutate(
    Ambang = factor(
      sprintf("|log2FC| >= %.2f\n(>= %.2fx lipat)", Threshold_abs_log2FC, 2^Threshold_abs_log2FC),
      levels = sprintf("|log2FC| >= %.2f\n(>= %.2fx lipat)", c(0.58, 1.00, 1.50), round(2^c(0.58, 1.00, 1.50), 2))
    ),
    Arah_Label = ifelse(Arah == "DEG_Up_in_Chitin",
                        "Naik pada Perlakuan Kitin (log2FC > 0)",
                        "Turun pada Perlakuan Kitin (log2FC < 0)")
  )

fig9 <- ggplot(sens_main, aes(x = Ambang, y = Jumlah, fill = Arah_Label)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Jumlah), position = position_dodge(width = 0.72),
            vjust = -0.4, fontface = "bold", size = 3.8) +
  scale_fill_manual(values = c("Naik pada Perlakuan Kitin (log2FC > 0)" = col_up,
                               "Turun pada Perlakuan Kitin (log2FC < 0)" = col_down)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Gambar 9. Jumlah DEG Naik vs Turun pada Berbagai Ambang Sensitivitas |log2FC| (padj < 0,05)",
    subtitle = "Kontras: Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam) | Estimator: DESeq2 + apeglm Shrinkage",
    x = "Ambang Batas |log2 Fold Change| (dengan FDR padj < 0,05)",
    y = "Jumlah Gen Ekspresi Diferensial (DEG)",
    fill = "Arah Regulasi (Perlakuan vs Kontrol)"
  ) +
  theme_bw(base_size = 11.5) +
  theme(legend.position = "top",
        axis.text.x = element_text(face = "bold", size = 10.5),
        plot.title = element_text(face = "bold", size = 12.5))
save_pub_plot(fig9, "Fig09_BarChart_DEG_Counts_and_Sensitivity", width = 9.5, height = 6)

cat("\n[SUKSES] Tahap 4 s.d. Tahap 8 selesai dijalankan.\n")
