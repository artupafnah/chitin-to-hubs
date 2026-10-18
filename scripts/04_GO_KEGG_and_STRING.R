#!/usr/bin/env Rscript
# ==============================================================================
# TAHAP 9: ANOTASI FUNGSIONAL (GO BP/MF/CC & KEGG) DAN FILE INPUT JEJARING STRING
# TAHAP 10: PENCATATAN VERSI PERANGKAT LUNAK & sessionInfo()
#
# Dataset   : GSE133053 / BioProject PRJNA549850 (Wan et al., 2019, Nature Plants)
# Organisme : Arabidopsis thaliana (Col-0, TAIR10, TaxID: 3702)
# Kontras   : Perlakuan Elisitor Kitin 1 jam (C6) vs Kontrol Air 1 jam (H2O)
# Latar     : Semua gen yang diuji dalam DESeq2 (universe = 19.979 gen)
# Ambang    : FDR (Benjamini-Hochberg p.adjust) < 0,05
# Seed      : 42
# ==============================================================================

set.seed(42)
options(stringsAsFactors = FALSE, timeout = 25)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(dplyr)
  library(tidyr)
})

cat("[TAHAP 9] Memulai analisis pengayaan fungsional GO (BP/MF/CC), KEGG, dan persiapan STRING...\n")

save_pub_plot <- function(plot_obj, base_name, width = 10, height = 7) {
  png_file <- file.path("figures", paste0(base_name, ".png"))
  pdf_file <- file.path("figures", paste0(base_name, ".pdf"))
  svg_file <- file.path("figures", paste0(base_name, ".svg"))
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_file, plot = plot_obj, width = width, height = height, bg = "white")
  tryCatch({
    ggsave(svg_file, plot = plot_obj, width = width, height = height, bg = "white")
  }, error = function(e) {
    svg(svg_file, width = width, height = height)
    print(plot_obj)
    dev.off()
  })
  cat(sprintf("[GAMBAR] Tersimpan: %s (.png 300 dpi, .pdf, .svg)\n", base_name))
}

# 1. Muat tabel hasil DESeq2 dari Tahap 7
res_df   <- read.csv("results/DESeq2_all_genes_complete_table.csv")
deg_all  <- read.csv("results/DEG_significant_all_logFC1.csv")
deg_up   <- read.csv("results/DEG_upregulated_in_Chitin_C6.csv")
deg_down <- read.csv("results/DEG_downregulated_in_Chitin_C6.csv")

universe_genes <- unique(res_df$gene_id)
up_genes       <- unique(deg_up$gene_id)
down_genes     <- unique(deg_down$gene_id)

cat(sprintf("[INFO] Latar belakang (universe): %d gen | DEG Naik: %d gen | DEG Turun: %d gen\n",
            length(universe_genes), length(up_genes), length(down_genes)))

# ==============================================================================
# 2. ANALISIS ENRICHMENT GENE ONTOLOGY (BP, MF, CC) TERPISAH UNTUK NAIK & TURUN
# ==============================================================================
run_go_all_ontologies <- function(gene_vec, direction_label, p_cut = 0.05, q_cut = 0.05) {
  out_list <- list()
  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("[GO] Menjalankan enrichGO (%s) untuk %s (%d gen)...\n", ont, direction_label, length(gene_vec)))
    ego <- enrichGO(
      gene          = gene_vec,
      universe      = universe_genes,
      OrgDb         = org.At.tair.db,
      keyType       = "TAIR",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cut,
      qvalueCutoff  = q_cut,
      readable      = FALSE
    )
    if (!is.null(ego)) {
      df <- ego@result
      df <- df[df$p.adjust <= p_cut, , drop = FALSE]
      if (nrow(df) > 0) {
        df$Ontology  <- ont
        df$Direction <- direction_label
        df$GeneRatio_num <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        out_list[[ont]] <- df
      }
    }
  }
  if (length(out_list) == 0) return(data.frame())
  bind_rows(out_list)
}

if (file.exists("results/GO_enrichment_UP_BP_MF_CC.csv") && file.size("results/GO_enrichment_UP_BP_MF_CC.csv") > 100) {
  cat("[GO] Memuat hasil enrichGO UP yang sudah dihitung dari results/GO_enrichment_UP_BP_MF_CC.csv...\n")
  go_up_df <- read.csv("results/GO_enrichment_UP_BP_MF_CC.csv")
} else {
  go_up_df <- run_go_all_ontologies(up_genes, "Naik pada Perlakuan Kitin (Up)", p_cut = 0.05, q_cut = 0.05)
  write.csv(go_up_df, "results/GO_enrichment_UP_BP_MF_CC.csv", row.names = FALSE)
}

# Jalankan untuk gen turun dan rekam baik FDR < 0.05 maupun p-value nominal < 0.05
down_list_fdr <- list()
down_list_nom <- list()
for (ont in c("BP", "MF", "CC")) {
  cat(sprintf("[GO] Menjalankan enrichGO (%s) untuk Turun pada Perlakuan Kitin (%d gen)...\n", ont, length(down_genes)))
  ego_d <- enrichGO(
    gene          = down_genes,
    universe      = universe_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1.0,
    qvalueCutoff  = 1.0,
    readable      = FALSE
  )
  if (!is.null(ego_d) && nrow(ego_d@result) > 0) {
    dfd <- ego_d@result
    dfd$Ontology  <- ont
    dfd$Direction <- "Turun pada Perlakuan Kitin (Down)"
    dfd$GeneRatio_num <- sapply(strsplit(dfd$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    dfd_fdr <- dfd[dfd$p.adjust < 0.05, , drop = FALSE]
    dfd_nom <- dfd[dfd$pvalue < 0.05, , drop = FALSE]
    if (nrow(dfd_fdr) > 0) down_list_fdr[[ont]] <- dfd_fdr
    if (nrow(dfd_nom) > 0) down_list_nom[[ont]] <- dfd_nom
  }
}
go_down_df  <- if (length(down_list_fdr) > 0) bind_rows(down_list_fdr) else data.frame(ID=character(), Description=character(), Ontology=character(), p.adjust=numeric())
go_down_nom <- if (length(down_list_nom) > 0) bind_rows(down_list_nom) else data.frame()
if (nrow(go_down_nom) > 0) go_down_nom <- go_down_nom[order(go_down_nom$pvalue), ]

write.csv(go_down_df,  "results/GO_enrichment_DOWN_BP_MF_CC.csv", row.names = FALSE)
write.csv(go_down_nom, "results/GO_enrichment_DOWN_nominal_diagnostic.csv", row.names = FALSE)

cat(sprintf("\n[GO] Jumlah term GO signifikan (FDR < 0,05): Naik = %d term, Turun = %d term (Nominal p<0,05 pada Turun = %d term)\n",
            nrow(go_up_df), nrow(go_down_df), nrow(go_down_nom)))

cat("\n=== Ringkasan Jumlah Term GO Signifikan (FDR < 0.05) pada Gen Naik per Kategori ===\n")
print(table(go_up_df$Ontology))

# Pilih Top term per ontologi untuk visualisasi Dotplot (Gambar 10)
top_go_up <- go_up_df %>%
  group_by(Ontology) %>%
  slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Term_Label = ifelse(nchar(Description) > 52,
                        paste0(ID, ": ", substr(Description, 1, 49), "..."),
                        paste0(ID, ": ", Description)),
    neg_log10_padj = -log10(p.adjust),
    Ontology_Full = case_when(
      Ontology == "BP" ~ "Biological Process (BP)",
      Ontology == "MF" ~ "Molecular Function (MF)",
      Ontology == "CC" ~ "Cellular Component (CC)"
    )
  ) %>%
  arrange(Ontology, GeneRatio_num) %>%
  mutate(Term_Label = factor(Term_Label, levels = unique(Term_Label)))

fig10 <- ggplot(top_go_up, aes(x = GeneRatio_num, y = Term_Label, size = Count, color = neg_log10_padj)) +
  geom_point(alpha = 0.9) +
  facet_grid(Ontology_Full ~ ., scales = "free_y", space = "free_y") +
  scale_color_viridis_c(option = "plasma", end = 0.9, name = "-log10(FDR)") +
  scale_size_continuous(range = c(3.5, 8.5), name = "Jumlah Gen DEG") +
  labs(
    title = "Gambar 10. Dotplot Pengayaan Fungsional Gene Ontology (BP & MF) pada DEG Naik (FDR < 0,05)",
    subtitle = sprintf("Kontras: Perlakuan Kitin C6 (1 jam) vs Kontrol H2O (1 jam) | Latar: %s Gen Teruji | Gen Turun: 0 term lolos FDR < 0,05",
                       format(length(universe_genes), big.mark = ".")),
    x = "Gene Ratio (Proporsi DEG terhadap Total Gen Beranotasi pada Grup)",
    y = "Term Gene Ontology (GO ID & Deskripsi)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#F0F4F8"),
    strip.text = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(color = "black", size = 9.5, face = "bold"),
    plot.title = element_text(face = "bold", size = 12.5)
  )
save_pub_plot(fig10, "Fig10_GO_Enrichment_Dotplot_Up_vs_Down", width = 11.5, height = 8)

# ==============================================================================
# 3. ANALISIS ENRICHMENT JALUR KEGG (TERPISAH NAIK & TURUN)
# ==============================================================================
cat("\n[KEGG] Menjalankan analisis pengayaan jalur metabolisme & sinyal KEGG...\n")

path_map <- AnnotationDbi::select(org.At.tair.db, keys = universe_genes, columns = "PATH", keytype = "TAIR")
path_map <- path_map[!is.na(path_map$PATH), c("PATH", "TAIR")]
path_map$PATH <- paste0("ath", path_map$PATH)

kegg_names <- c(
  "ath04626" = "Plant-pathogen interaction",
  "ath04016" = "MAPK signaling pathway - plant",
  "ath04075" = "Plant hormone signal transduction",
  "ath00940" = "Phenylpropanoid biosynthesis",
  "ath00900" = "Terpenoid backbone biosynthesis",
  "ath00270" = "Cysteine and methionine metabolism (Ethylene)",
  "ath00480" = "Glutathione metabolism",
  "ath00592" = "alpha-Linolenic acid metabolism (Jasmonate)",
  "ath00966" = "Glucosinolate biosynthesis",
  "ath00195" = "Photosynthesis",
  "ath00196" = "Photosynthesis - antenna proteins",
  "ath00500" = "Starch and sucrose metabolism",
  "ath00520" = "Amino sugar and nucleotide sugar metabolism",
  "ath04141" = "Protein processing in endoplasmic reticulum",
  "ath00910" = "Nitrogen metabolism",
  "ath00904" = "Diterpenoid biosynthesis",
  "ath00906" = "Carotenoid biosynthesis (ABA precursor)",
  "ath00360" = "Phenylalanine metabolism",
  "ath00400" = "Phenylalanine, tyrosine and tryptophan biosynthesis"
)
u_paths <- unique(path_map$PATH)
term2name <- data.frame(
  PATH = u_paths,
  Name = ifelse(u_paths %in% names(kegg_names), kegg_names[u_paths], paste0("KEGG Pathway ", u_paths))
)

run_kegg_analysis <- function(gene_vec, direction_label) {
  # Coba enrichKEGG daring terlebih dahulu
  res_online <- tryCatch({
    ek <- enrichKEGG(
      gene          = gene_vec,
      universe      = universe_genes,
      organism      = "ath",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1.0,
      qvalueCutoff  = 1.0
    )
    if (!is.null(ek) && nrow(ek@result) > 0) {
      df <- ek@result
      df$Description <- sub(" - Arabidopsis thaliana \\(thale cress\\)$", "", df$Description)
      df$Direction <- direction_label
      df$Source <- "KEGG REST API (ath)"
      return(df)
    } else NULL
  }, error = function(e) {
    cat("[INFO] Beralih ke anotasi lokal org.At.tairPATH:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(res_online)) return(res_online)

  ek_loc <- enricher(
    gene          = gene_vec,
    universe      = universe_genes,
    TERM2GENE     = path_map,
    TERM2NAME     = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1.0,
    qvalueCutoff  = 1.0
  )
  if (!is.null(ek_loc) && nrow(ek_loc@result) > 0) {
    df <- ek_loc@result
    df$Direction <- direction_label
    df$Source <- "org.At.tair.db Local KEGG"
    return(df)
  }
  data.frame()
}

kegg_up_full   <- run_kegg_analysis(up_genes,   "Naik pada Perlakuan Kitin (Up)")
kegg_down_full <- run_kegg_analysis(down_genes, "Turun pada Perlakuan Kitin (Down)")

kegg_up_sig   <- kegg_up_full[kegg_up_full$p.adjust < 0.05, , drop = FALSE]
kegg_down_sig <- kegg_down_full[kegg_down_full$p.adjust < 0.05, , drop = FALSE]
kegg_down_nom <- kegg_down_full[kegg_down_full$pvalue < 0.05, , drop = FALSE]

kegg_sig_all <- bind_rows(kegg_up_sig, kegg_down_sig)
write.csv(kegg_sig_all,  "results/KEGG_enrichment_UP_and_DOWN.csv", row.names = FALSE)
write.csv(kegg_down_nom, "results/KEGG_enrichment_DOWN_nominal_diagnostic.csv", row.names = FALSE)

cat(sprintf("\n[KEGG] Jalur KEGG signifikan (FDR < 0,05): Naik = %d jalur, Turun = %d jalur (Nominal p<0,05 pada Turun = %d jalur)\n",
            nrow(kegg_up_sig), nrow(kegg_down_sig), nrow(kegg_down_nom)))
if (nrow(kegg_sig_all) > 0) {
  print(kegg_sig_all[, c("ID", "Description", "Direction", "GeneRatio", "pvalue", "p.adjust", "Count")])
}
if (nrow(kegg_down_nom) > 0) {
  cat("\n=== Jalur KEGG pada Gen Turun (Nominal p < 0.05) ===\n")
  print(kegg_down_nom[, c("ID", "Description", "Direction", "GeneRatio", "pvalue", "p.adjust", "Count")])
}

# Buat plot gabungan KEGG (Gambar 11): Jalur signifikan pada Gen Naik (FDR < 0.05) + Jalur teratas Gen Turun
plot_kegg_list <- list()
if (nrow(kegg_up_sig) > 0) {
  ku <- head(kegg_up_sig[order(kegg_up_sig$p.adjust), ], 10)
  ku$Group_Panel <- "DEG Naik pada Kitin (FDR < 0,05)"
  plot_kegg_list[[1]] <- ku
} else if (nrow(kegg_up_full) > 0) {
  ku <- head(kegg_up_full[order(kegg_up_full$pvalue), ], 8)
  ku$Group_Panel <- "DEG Naik pada Kitin (Top Jalur)"
  plot_kegg_list[[1]] <- ku
}
if (nrow(kegg_down_sig) > 0) {
  kd <- head(kegg_down_sig[order(kegg_down_sig$p.adjust), ], 5)
  kd$Group_Panel <- "DEG Turun pada Kitin (FDR < 0,05)"
  plot_kegg_list[[2]] <- kd
} else if (nrow(kegg_down_nom) > 0) {
  kd <- head(kegg_down_nom[order(kegg_down_nom$pvalue), ], 5)
  kd$Group_Panel <- "DEG Turun pada Kitin (Nominal p < 0,05)"
  plot_kegg_list[[2]] <- kd
}

plot_kegg_df <- bind_rows(plot_kegg_list)
if (nrow(plot_kegg_df) > 0) {
  plot_kegg_df$neg_log10_padj <- -log10(plot_kegg_df$p.adjust)
  plot_kegg_df$Path_Label <- paste0(plot_kegg_df$ID, ": ", plot_kegg_df$Description)
  plot_kegg_df <- plot_kegg_df %>%
    arrange(Group_Panel, Count) %>%
    mutate(Path_Label = factor(Path_Label, levels = unique(Path_Label)))

  fig11 <- ggplot(plot_kegg_df, aes(x = Count, y = Path_Label, fill = neg_log10_padj)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = sprintf("FDR=%.2e (p=%.2e, n=%d)", p.adjust, pvalue, Count)),
              hjust = -0.04, size = 3.0, fontface = "bold") +
    facet_grid(Group_Panel ~ ., scales = "free_y", space = "free_y") +
    scale_fill_viridis_c(option = "magma", begin = 0.25, end = 0.9, name = "-log10(FDR)") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
    labs(
      title = "Gambar 11. Barplot Pengayaan Jalur KEGG (KEGG Pathway Enrichment): Perlakuan Kitin C6 vs Kontrol H2O",
      subtitle = "Latar: 19.979 Gen Teruji (Arabidopsis thaliana, ath) | Kontras: Perlakuan Kitin C6 (1 jam) vs Kontrol H2O (1 jam)",
      x = "Jumlah Gen DEG dalam Jalur KEGG",
      y = "Jalur Metabolisme / Transduksi Sinyal KEGG"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#F0F4F8"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.y = element_text(color = "black", face = "bold", size = 9.5),
      plot.title = element_text(face = "bold", size = 12.5)
    )
  save_pub_plot(fig11, "Fig11_KEGG_Pathway_Enrichment_Plot", width = 12.5, height = 6.8)
}

# ==============================================================================
# 4. PERSIAPAN BERKAS INPUT UNTUK JEJARING PPI STRING (TAHAP 9)
# ==============================================================================
cat("\n[STRING] Menyiapkan berkas input untuk analisis jejaring interaksi protein STRING...\n")

unique_tair_deg <- unique(deg_all$gene_id)
writeLines(unique_tair_deg, "results/STRING_DEG_unique_TAIR_IDs.txt")

unique_sym_deg <- unique(deg_all$symbol)
writeLines(unique_sym_deg, "results/STRING_DEG_unique_GeneSymbols.txt")

string_val_df <- deg_all %>%
  dplyr::select(gene_id, log2FC) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(desc(log2FC))

write.table(string_val_df, "results/STRING_DEG_with_log2FC.tsv",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

string_cyto_df <- deg_all %>%
  dplyr::select(gene_id, symbol, log2FC, padj, regulation, description) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(desc(abs(log2FC)))

write.table(string_cyto_df, "results/STRING_Cytoscape_node_attributes_DEG.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat(sprintf("[STRING] Berkas siap unggah tersimpan di results/ (%d gen unik DEG):\n", length(unique_tair_deg)))
cat("  1. results/STRING_DEG_unique_TAIR_IDs.txt (1 ID TAIR per baris untuk 'Multiple proteins')\n")
cat("  2. results/STRING_DEG_unique_GeneSymbols.txt (1 Simbol Gen per baris)\n")
cat("  3. results/STRING_DEG_with_log2FC.tsv (Pasangan ID TAIR dan nilai log2FC)\n")
cat("  4. results/STRING_Cytoscape_node_attributes_DEG.tsv (Tabel atribut lengkap untuk Cytoscape)\n")

# ==============================================================================
# 5. SIMPAN sessionInfo() DAN VERSI PERANGKAT LUNAK KE logs/
# ==============================================================================
sink("logs/sessionInfo_and_software_versions.txt")
cat("==============================================================================\n")
cat("INFORMASI SESI R DAN VERSI PERANGKAT LUNAK (REPRODUKSIBILITAS PIPELINE)\n")
cat("Tanggal Eksekusi :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Random Seed      : 42\n")
cat("Dataset          : NCBI GEO GSE133053 / BioProject PRJNA549850 / SRA SRP201971\n")
cat("Organisme        : Arabidopsis thaliana (Col-0, TAIR10)\n")
cat("Kontras Utama    : Perlakuan Kitin C6 (1 jam) vs Kontrol Air H2O (1 jam)\n")
cat("==============================================================================\n\n")
print(sessionInfo())
cat("\n--- Versi Paket Utama R ---\n")
main_pkgs <- c("DESeq2", "apeglm", "clusterProfiler", "org.At.tair.db", "AnnotationDbi",
               "ggplot2", "ggrepel", "pheatmap", "viridis", "patchwork", "dplyr", "tidyr")
for (p in main_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("%-18s : %s\n", p, as.character(packageVersion(p))))
  }
}
sink()

cat("\n[SUKSES] Tahap 9 dan pencatatan sessionInfo() selesai dijalankan.\n")
