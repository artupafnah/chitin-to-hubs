# chitin-to-hubs

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R: 4.5.2](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)
[![Python: 3.11](https://img.shields.io/badge/Python-3.11-brightgreen.svg)](https://www.python.org/)
[![Bioconductor: 3.22](https://img.shields.io/badge/Bioconductor-3.22-green.svg)](https://bioconductor.org/)

> **From chitin elicitation to PPI hubs: RNA-seq differential expression, enrichment, and STRING network analysis of GSE133053 (chitin vs water, 1 h).**

Repositori ini memuat seluruh sumber daya data kuantitatif, skrip analisis terotomatisasi (*pipeline*), berkas log eksekusi, serta visualisasi beresolusi publikasi (PNG 300 dpi, PDF vektor, dan SVG) untuk rekonstruksi komprehensif respons dini kekebalan tanaman (*Pattern-Triggered Immunity* / PTI) pada *Arabidopsis thaliana* pasca-paparan elisitor kitin heksamer ($\text{C}_6$).

---

## Daftar Isi
- [Tujuan Penelitian](#tujuan-penelitian)
- [Dataset & Desain Eksperimen](#dataset--desain-eksperimen)
- [Arsitektur Pipeline Analisis](#arsitektur-pipeline-analisis)
- [Struktur Direktori Repositori](#struktur-direktori-repositori)
- [Panduan Menjalankan Ulang (Reproducibility)](#panduan-menjalankan-ulang-reproducibility)
- [Ringkasan Temuan Utama](#ringkasan-temuan-utama)
- [Kandidat Hub Prioritas (Top 10)](#kandidat-hub-prioritas-top-10)
- [Keterbatasan Metodologis](#keterbatasan-metodologis)
- [Sumber Data & Sitasi](#sumber-data--sitasi)
- [Lisensi & Penulis](#lisensi--penulis)

---

## Tujuan Penelitian
1. Mengevaluasi profil kendali mutu (*quality control* / QC) kedalaman pustaka sekuensing RNA-seq dan mendeteksi anomali teknis/sampel *outlier* pada dataset publik NCBI GEO **GSE133053**.
2. Mengidentifikasi gen-gen yang berdiferensiasi ekspresi secara signifikan (*Differentially Expressed Genes* / DEG) antara perlakuan elisitor kitin heksamer ($\text{C}_6$, 1 jam) dan kontrol air ($\text{H}_2\text{O}$, 1 jam) menggunakan model GLM Binomial Negatif `DESeq2` dan penyusutan Bayesian empiris berprior Cauchy `apeglm`.
3. Memetakan pengayaan fungsional proses biologis (GO:BP), fungsi molekuler (GO:MF), komponen seluler (GO:CC), dan jalur metabolisme/pensinyalan KEGG secara terpisah antara gen naik (*upregulated*) dan gen turun (*downregulated*).
4. Merekonstruksi jaringan interaksi protein (*Protein-Protein Interaction* / PPI) berbasis basis data STRING v12.0, mengevaluasi topologi bebas-skala (*scale-free*), mempartisi modul komunitas fungsional dengan algoritma Louvain, serta memprioritaskan sepuluh kandidat gen *hub* utama menggunakan Skor Prioritas Komposit empat metrik sentralitas (*Degree*, *Betweenness*, *Closeness*, *Eigenvector*).

---

## Dataset & Desain Eksperimen
- **Aksesi Repositori Publik:** NCBI GEO [GSE133053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133053) | BioProject [PRJNA549850](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA549850) | SRA SRP201971.
- **Organisme:** *Arabidopsis thaliana* ekotipe Columbia-0 (Col-0, NCBI Taxonomy ID: `3702`).
- **Publikasi Sumber:** Wan *et al.* (2019), *Nature Plants* 5(2):176–187 ([DOI: 10.1038/s41477-018-0339-8](https://doi.org/10.1038/s41477-018-0339-8)).
- **Desain Eksperimen Subset:**
  - Jaringan kecambah usia 12 hari yang diinkubasi semalam (*overnight resting*) dalam air steril.
  - Perlakuan: Elisitor kitin heksamer terpurifikasi ($\text{C}_6$, konsentrasi 100 µM) selama 1 jam (*n* = 4 replikat biologis: `C6_1h_Rep1` s.d. `Rep4`).
  - Kontrol: Air steril bebas elisitor ($\text{H}_2\text{O}$) selama 1 jam (*n* = 4 replikat biologis: `H2O_1h_Rep1` s.d. `Rep4`).
  - Platform Sekuensing: Illumina HiSeq 2000, *Single-End* (SE) 101 bp, pustaka cDNA *unstranded poly-A selection*.
- **Arah Kontras Analisis:**
  $$\text{Kontras} = \text{Perlakuan Kitin } \text{C}_6\ (1\text{ jam}) \quad\text{vs}\quad \text{Kontrol Air } \text{H}_2\text{O}\ (1\text{ jam})$$
  - $\log_2\text{FC} > 0$: Gen mengalami **peningkatan ekspresi (Naik / *Upregulated*)** akibat paparan kitin.
  - $\log_2\text{FC} < 0$: Gen mengalami **penurunan ekspresi (Turun / *Downregulated*)** akibat paparan kitin.

---

## Arsitektur Pipeline Analisis

```
[NCBI SRA / GEO GSE133053]
       │
       ▼ (scripts/01_download_metadata_and_counts.py)
[Metadata & Arsip GSE133053_RAW.tar]
       │
       ▼ (scripts/02_extract_and_build_count_matrix.py)
[Ekstraksi 52 GSM TSV -> Matriks Count 31.591 Lokus TAIR10]
       │  Audit QC Pustaka: C6_1h_Rep1 gagal (0,64% mapping) -> Dikeluarkan (qc/)
       ▼
[Matriks Bersih 7 Sampel Valid (4 Kontrol vs 3 Perlakuan)]
       │
       ▼ (scripts/03_DESeq2_DEG_and_Visualizations.R)
[Penyaringan Awal: 19.979 Gen Aktif (count >= 10 pada >= 3 sampel)]
       │  Normalisasi: Median-of-Ratios DESeq2 & VST
       │  Pemodelan: GLM Binomial Negatif & Wald Test
       │  Penyusutan Efek: apeglm (Prior Cauchy)
       ▼
[717 DEG Signifikan: padj < 0,05 & |log2FC| >= 1,0]
(616 Gen Naik / 85,91%  vs  101 Gen Turun / 14,09%)
       │
       ├─────────────────────────────────────────┐
       ▼ (scripts/04 & 07)                       ▼ (scripts/05_prepare_deg_lists_and_string_mapping.py)
[Analisis Pengayaan ORA]                 [Pemetaan Proteom STRING v12.0 API]
ClusterProfiler (GO:BP/MF/CC & KEGG)     690 Protein Terpetakan (96,23%)
All DEG, Up DEG, Down DEG, Top 10 Hub     27 Lokus lncRNA/Pseudogen Tak Terpetakan
       │                                         │
       │                                         ▼ (scripts/06_string_ppi_network_centrality_and_prioritization.py)
       │                                 [Rekonstruksi Jaringan PPI STRING]
       │                                 GCC: 488 Simpul, 1.513 Sisi (Skor >= 0,400)
       │                                         │
       │                                         ├─ Analisis Topologi Bebas-Skala (gamma = 1,44; R^2 = 0,89)
       │                                         ├─ Partisi Modul Komunitas Louvain (Q = 0,5603; 16 Modul)
       │                                         ├─ Sentralitas: Degree, Betweenness, Closeness, Eigenvector
       │                                         └─ Standarisasi Z-score (STDEV.S) -> Skor Prioritas Komposit
       ▼                                         ▼
[Visualisasi Publikasi 300 dpi (figures/)]  [10 Kandidat Hub Imunitas Prioritas (results/)]
```

---

## Struktur Direktori Repositori

```
chitin-to-hubs/
├── data/
│   ├── raw_counts_gsm/                 # 52 berkas .txt.gz kuantifikasi eXpress per GSM
│   ├── GSE133053_RAW.tar               # Arsip asli matriks count resmi GEO (7,1 MB)
│   ├── GSE133053_gsm_brief.soft        # Ekstraksi metadata SOFT resmi GEO
│   ├── metadata.csv                    # Metadata 8 sampel titik waktu 1 jam
│   ├── metadata_valid7samples.csv      # Metadata 7 sampel lolos kendali mutu
│   ├── metadata_GSE133053_all52samples.csv # Metadata lengkap 52 sampel seri GSE133053
│   └── SraRunTable_PRJNA549850.csv     # Metadata run SRA dari Run Selector NCBI
├── qc/
│   └── qc_library_and_mapping_summary.csv # Tabel evaluasi penugasan read & status QC
├── scripts/
│   ├── 00_upstream_fastq_to_counts_reference.sh # Skrip referensi FastQC/fastp/HISAT2/featureCounts
│   ├── 01_download_metadata_and_counts.py       # Pengunduhan metadata NCBI & verifikasi berkas
│   ├── 02_extract_and_build_count_matrix.py     # Ekstraksi matriks count & audit pustaka
│   ├── 03_DESeq2_DEG_and_Visualizations.R       # Normalisasi, DESeq2, apeglm, PCA, plot DEG
│   ├── 04_GO_KEGG_and_STRING.R                  # Pengayaan GO/KEGG clusterProfiler awal
│   ├── 05_prepare_deg_lists_and_string_mapping.py # Kurasi daftar DEG & pemetaan API STRING
│   ├── 06_string_ppi_network_centrality_and_prioritization.py # PPI, Louvain, sentralitas & skor komposit
│   └── 07_enrichment_all_up_down_and_top10.R    # Pengayaan komparatif All, Up, Down, Top 10
├── results/
│   ├── DESeq2_all_genes_complete_table.csv      # Tabel lengkap hasil pengujian 19.979 gen aktif
│   ├── DEG_significant_all_logFC1.csv           # 717 DEG utama (|log2FC| >= 1,0; padj < 0,05)
│   ├── DEG_upregulated_in_Chitin_C6.csv         # 616 DEG Naik
│   ├── DEG_downregulated_in_Chitin_C6.csv       # 101 DEG Turun
│   ├── DEG_sensitivity_summary_table.csv        # Uji sensitivitas ambang log2FC & sampel
│   ├── normalized_counts_DESeq2.csv             # Matriks count ternormalisasi median-of-ratios
│   ├── vst_transformed_matrix.csv               # Matriks ekspresi tertransformasi VST
│   ├── PCA_coordinates_valid7samples.csv        # Koordinat PC1 s.d. PC4 pada 7 sampel valid
│   ├── sample_correlation_matrix_VST.csv        # Matriks korelasi Pearson antarsampel
│   ├── sample_euclidean_distance_VST.csv        # Matriks jarak Euclidean antarsampel
│   ├── enrichment_ALL_DEG_GO_KEGG.csv           # Tabel ORA lengkap 717 DEG
│   ├── enrichment_UP_DEG_GO_KEGG.csv            # Tabel ORA 616 DEG Naik
│   ├── enrichment_DOWN_DEG_GO_KEGG.csv          # Tabel ORA 101 DEG Turun
│   ├── enrichment_TOP10_hubs_GO_KEGG.csv        # Tabel ORA 10 gen hub prioritas
│   ├── enrichment_comparison_Top10_vs_AllDEG.csv# Perbandingan ORA Top 10 Hub vs All DEG
│   ├── STRING_interactions_score0400_medium.csv # 1.522 sisi interaksi PPI (skor >= 0,400)
│   ├── STRING_interactions_score0700_high.csv   # 284 sisi interaksi PPI (skor >= 0,700)
│   ├── PPI_centrality_GCC_score0400_complete.csv# Sentralitas lengkap 488 simpul GCC
│   ├── PPI_Louvain_communities_summary.csv      # Rangkuman 16 modul komunitas Louvain
│   ├── top_candidates.csv                       # Sepuluh kandidat hub prioritas teratas
│   ├── top_candidates.xlsx                      # Buku kerja Excel multi-lembar terformat
│   └── PPI_top10_stability_comparison.csv       # Uji kestabilan peringkat melintasi 5 skema
├── figures/                                     # 18 set visualisasi (PNG 300 dpi, PDF, SVG)
│   ├── Fig01_QC_Library_and_Mapping_Rate.*
│   ├── Fig02_Boxplot_Before_vs_After_Normalization.*
│   ├── Fig03_PCA_Plot_VST_Top500.*
│   ├── Fig04_Sample_Correlation_Distance_Heatmap.*
│   ├── Fig05_Volcano_Plot_C6_vs_H2O.*
│   ├── Fig06_MA_Plot_Before_and_After_apeglm.*
│   ├── Fig07_Heatmap_Top50_DEG_Zscore.*
│   ├── Fig08_Pvalue_Distribution_Histogram.*
│   ├── Fig09_BarChart_DEG_Counts_and_Sensitivity.*
│   ├── Fig10_GO_Enrichment_Dotplot_Up_vs_Down.*
│   ├── Fig11_KEGG_Pathway_Enrichment_Plot.*
│   ├── Fig_PPI_01_Enrichment_All_DEG_GO_KEGG.*
│   ├── Fig_PPI_02_Enrichment_Up_vs_Down_Comparison.*
│   ├── Fig_PPI_03_STRING_PPI_Network_Overview.*
│   ├── Fig_PPI_04_Louvain_Community_and_Hub_Subnetwork.*
│   ├── Fig_PPI_05_Degree_Distribution_and_Topology.*
│   ├── Fig_PPI_06_Top10_Centrality_and_Composite_Score.*
│   └── Fig_PPI_07_Top10_vs_AllDEG_Enrichment_Comparison.*
├── logs/                                        # Catatan log eksekusi dan versi perangkat lunak
│   ├── 02_extract_and_build_count_matrix.log
│   ├── 03_DESeq2_DEG_and_Visualizations.log
│   ├── 04_GO_KEGG_and_STRING.log
│   ├── 05_prepare_deg_lists_and_string_mapping.log
│   ├── 06_string_ppi_network_centrality_and_prioritization.log
│   ├── 07_enrichment_all_up_down_and_top10.log
│   └── sessionInfo_and_software_versions.txt
├── requirements.txt                             # Dependensi pustaka Python
├── install_packages.R                           # Skrip instalasi dependensi R & Bioconductor
├── .gitignore
└── LICENSE                                      # Lisensi MIT
```

---

## Panduan Menjalankan Ulang (Reproducibility)

### 1. Kebutuhan Lingkungan Perangkat Lunak
- **Sistem Operasi:** Windows 10/11 x64 atau Linux x86_64
- **Python:** Versi ≥ 3.10 (teruji pada Python 3.11.9)
- **R:** Versi ≥ 4.3 (teruji pada R 4.5.2 / Bioconductor 3.22)
- **Koneksi Internet:** Diperlukan untuk instalasi paket dan kueri REST API basis data STRING v12.0

### 2. Instalasi Dependensi
```bash
# Instal dependensi Python
pip install -r requirements.txt

# Instal dependensi R & Bioconductor
Rscript install_packages.R
```

### 3. Eksekusi Skrip Berurutan
Seluruh skrip analitik telah dikonfigurasi dengan *fixed random seed* (`seed = 42`) dan menggunakan *relative paths* terhadap *root* repositori:

```bash
# Tahap 1: Ekstraksi matriks count & audit kendali mutu sampel
python scripts/02_extract_and_build_count_matrix.py

# Tahap 2: Normalisasi DESeq2, penyusutan apeglm, dan visualisasi DEG
Rscript scripts/03_DESeq2_DEG_and_Visualizations.R

# Tahap 3: Persiapan daftar lokus unik & pemetaan ke proteom STRING v12.0
python scripts/05_prepare_deg_lists_and_string_mapping.py

# Tahap 4: Pemodelan jaringan PPI STRING, modularitas Louvain, sentralitas & prioritisasi hub
python scripts/06_string_ppi_network_centrality_and_prioritization.py

# Tahap 5: Analisis pengayaan fungsional komparatif (All DEG, Up, Down, Top 10 Hub)
Rscript scripts/07_enrichment_all_up_down_and_top10.R
```

---

## Ringkasan Temuan Utama

1. **Audit Mutu Pustaka Sekuensing (Tabel QC):**
   - Enam sampel menunjukkan kualitas pustaka sangat baik (10,88–19,99 juta *reads* SRA, penugasan unik 87,23%–93,19%).
   - Sampel `C6_1h_Rep1` (`GSM3899022`) teridentifikasi mengalami kegagalan teknis sekuensing parah (hanya 4.816 *reads* terpetakan atau 0,64%) dan dikeluarkan dari pemodelan utama (*FAIL_EXCLUDE*).
   - Sampel `C6_1h_Rep3` (`GSM3899036`) memiliki persentase penugasan 12,29% namun memuat 2.375.723 *reads* mRNA valid pada 15.935 gen aktif; normalisasi *median-of-ratios* ($s_j = 0,2124$) berhasil mengangkat ukuran efektifnya menjadi 11,18 juta *reads*, sejajar dengan sampel lainnya.
2. **Pemisahan Transkriptomik Global (PCA & Korelasi):**
   - Eksklusi `C6_1h_Rep1` mengembalikan orientasi sumbu biologis pada PCA: PC1 menjelaskan **72,8% varians** dan memisahkan secara tegas kelompok perlakuan kitin ($\text{C}_6$) terhadap kontrol air ($\text{H}_2\text{O}$).
   - Korelasi Pearson internal kelompok perlakuan mencapai $r = 0,936$–$0,969$, dan kontrol mencapai $r = 0,974$–$0,983$.
3. **Asimetri Ekstrem Respons Imun Kitin:**
   - Dari 19.979 gen aktif, teridentifikasi **717 DEG** ($\text{FDR} < 0,05$ dan $|\log_2\text{FC}_{\text{shrunk}}| \ge 1,0$).
   - Sebanyak **616 gen (85,91%) mengalami peningkatan ekspresi (Naik)**, dipimpin oleh enzim katabolisme asam absisat *CYP707A3* ($\log_2\text{FC} = +4,64; \text{padj} = 1,26 \times 10^{-63}$), reseptor kinase *CRK10* ($+3,54$) dan *CRK41* ($+3,58$), kaperon RE *AtCRT3* ($+2,02$), serta protein lignifikasi dinding sel *DIR7-like* ($+5,55$).
   - Hanya **101 gen (14,09%) mengalami penurunan ekspresi (Turun)**, merefleksikan penghentian asimilasi hara dan pemanjangan sel (*growth–defense trade-off*): transporter amonium *AMT1;3* ($\log_2\text{FC} = -2,04$), pengangkut nitrat *AtNPF7.3* ($-1,96$), kanal anion *SLAH4* ($-3,26$), dan peroksidase pemanjangan dinding sel *PER21* ($-1,70$).
4. **Pengayaan Fungsional Jalur Kanonik Imunitas:**
   - Kelompok DEG Naik memperkaya secara sangat kuat jalur KEGG **`ath04626` (*Plant-pathogen interaction*, 35 gen, $\text{FDR} = 8,33 \times 10^{-19}$)** dan **`ath04016` (*MAPK signaling pathway - plant*, 11 gen, $\text{FDR} = 0,0321$)**, serta istilah GO:BP respon imun bawaan (`GO:0045087`, 41 gen) dan transduksi sinyal kinase reseptor (`GO:0002221`, Fold Enrichment = 11,50×).
5. **Topologi Bebas-Skala dan Partisi Komunitas Louvain:**
   - Jaringan PPI STRING ($\text{skor} \ge 0,400$) pada komponen terhubung terbesar (GCC, 488 simpul, 1.513 interaksi) mengikuti distribusi hukum pangkat bebas-skala (*scale-free*, $\gamma = 1,44; R^2 = 0,89$).
   - Algoritma Louvain mempartisi jaringan ke dalam 16 modul modularitas tinggi ($Q = 0,5603$): Modul M1 (persepsi reseptor *AtLYK5*, kaskade *ATMPK11*, sinyal *ATEDS1*), Modul M2 (sensor kalmodulin *ATCAL5* & faktor transkripsi bZIP/WRKY), Modul M3 (regulasi umpan balik negatif E3 ubiquitin ligase *PUB24* & *CML40*), dan Modul M4 (eksekutor ledakan oksidatif *ATRBOHD* & kinase kalsium *CPK28*).

---

## Kandidat Hub Prioritas (Top 10)

Berdasarkan integrasi Skor Prioritas Komposit (*Composite Priority Score* = rata-rata empat *Z-score* sentralitas pada GCC), kesepuluh kandidat teratas seluruhnya merupakan gen yang mengalami induksi ekspresi signifikan ($\log_2\text{FC} > 0$) dan merekonstruksi kaskade transduksi sinyal kitin secara utuh:

| Rank | Lokus TAIR10 | Simbol Gen | $\log_2\text{FC}$ (`apeglm`) | Nilai `padj` | Modul Louvain | Derajat ($k$) | $Z_{\text{Deg}}$ | $Z_{\text{Bet}}$ | $Z_{\text{Clo}}$ | $Z_{\text{Eig}}$ | Skor Komposit | Peran Biologis Utama & Referensi Literatur |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **#1** | `AT5G47910` | **ATRBOHD** | +1,22 | $3,99 \times 10^{-3}$ | M4 | 52 | +6,30 | +9,55 | +2,83 | +6,05 | **+6,1812** | Enzim NADPH oxidase membran penghasil ledakan ROS apoplas saat elisitasi (Kadota *et al.* 2014). |
| **#2** | `AT5G66210` | **CPK28** | +2,42 | $5,87 \times 10^{-5}$ | M4 | 45 | +5,34 | +3,99 | +2,30 | +5,66 | **+4,3228** | Calcium-Dependent Protein Kinase penyangga imun pengatur pergantian protein BIK1 (Monaghan *et al.* 2014). |
| **#3** | `AT3G48090` | **ATEDS1** | +1,05 | $2,34 \times 10^{-2}$ | M1 | 33 | +3,69 | +5,13 | +2,23 | +3,81 | **+3,7138** | Node konvergensi esensial penghubung sinyal reseptor permukaan dengan imunitas intraseluler (Pruitt *et al.* 2021). |
| **#4** | `AT1G01560` | **ATMPK11** | +2,57 | $1,61 \times 10^{-2}$ | M1 | 37 | +4,24 | +2,75 | +2,20 | +5,34 | **+3,6302** | Kinase kaskade MAPK terinduksi-transkripsi pengirim sinyal ke faktor transkripsi nukleus (Bethke *et al.* 2012). |
| **#5** | `AT2G41110` | **ATCAL5** (*CAM2*) | +1,46 | $1,30 \times 10^{-2}$ | M2 | 33 | +3,69 | +5,29 | +2,37 | +3,15 | **+3,6243** | Isoform kalmodulin sensor kalsium sitosol pengaktivasi faktor transkripsi respon patogen (Cheval *et al.* 2013). |
| **#6** | `AT2G33580` | **AtLYK5** | +1,95 | $8,45 \times 10^{-4}$ | M1 | 37 | +4,24 | +2,89 | +2,09 | +4,91 | **+3,5308** | Reseptor kinase pengikat kitin berafinitas tinggi primer di membran plasma (Cao *et al.* 2014). |
| **#7** | `AT3G01830` | **CML40** | +2,24 | $1,69 \times 10^{-2}$ | M3 | 34 | +3,82 | +3,19 | +2,16 | +4,28 | **+3,3617** | Calmodulin-Like protein sensor $\text{Ca}^{2+}$ sitosol; kandidat regulator baru respon infeksi. |
| **#8** | `AT3G11840` | **PUB24** | +1,62 | $4,76 \times 10^{-2}$ | M3 | 32 | +3,55 | +4,00 | +2,11 | +3,55 | **+3,3014** | Enzim E3 ubiquitin ligase peredam respon imun melalui degradasi reseptor (Trujillo *et al.* 2008). |
| **#9** | `AT2G41100` | **ATCAL4** (*TCH3*) | +1,77 | $2,63 \times 10^{-2}$ | M1 | 30 | +3,27 | +2,67 | +2,01 | +3,16 | **+2,7783** | Kalmodulin sensor kalsium berulang internal terkait penguatan respon stres mekanik & biotik. |
| **#10** | `AT1G68880` | **AtbZIP** (*bZIP8*) | +1,99 | $1,49 \times 10^{-4}$ | M2 | 33 | +3,69 | +3,69 | +1,89 | +1,75 | **+2,7522** | Faktor transkripsi keluarga bZIP pengatur ekspresi gen hilir pertahanan sel. |

---

## Keterbatasan Metodologis
1. **Ukuran Replikat Efektif ($n = 4$ Kontrol vs $n = 3$ Perlakuan):** Eksklusi sampel gagal `C6_1h_Rep1` (0,64% *mapping rate*) mengurangi derajat bebas kelompok perlakuan kitin. Sampel `C6_1h_Rep3` memiliki kedalaman mentah lebih dangkal (2,38 juta *reads* terpetakan) meskipun berhasil dinormalisasi secara stabil ($r = 0,936$–$0,954$).
2. **Potret Satu Titik Waktu (1 Jam):** Analisis ini membedah fase puncak induksi transkripsi awal kitin dan belum memetakan kinetika peredaman jangka panjang pada 6 jam atau 24 jam.
3. **Kesenjangan Transkriptomik vs Fosfoproteomik:** Protein penting yang diaktivasi melalui fosforilasi cepat murni tanpa perubahan transkripsi mRNA $\ge 2$ kali lipat pada 1 jam (seperti ko-reseptor kinase *CERK1* `AT3G21630` atau kinase *MPK3*/*MPK6*) tidak muncul dalam daftar DEG.
4. **Bias Literatur Basis Data STRING:** STRING cenderung memberikan skor interaksi lebih tinggi pada protein yang telah diteliti intensif. Terdapat 27 lokus DEG (3,77%, sebagian besar lncRNA atau pseudogen) yang belum memiliki entri protein pada STRING v12.0.
5. **Karakter Hipotesis:** Peringkat sentralitas merupakan prioritisasi komputasi yang memerlukan konfirmasi eksperimental laboratorium basah (*qRT-PCR*, mutan T-DNA, uji interaksi *BiFC*/*Co-IP*).

---

## Sumber Data & Sitasi

### Repositori Data Primer
- NCBI GEO Series: [GSE133053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133053)
- NCBI BioProject: [PRJNA549850](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA549850)

### Sitasi Terkait
1. **Wan WL, Zhang L, Pruitt R, Zaidem M, Brugman R, Ma X, Krol E, Perraki A, Kilian J, Grossmann G, *et al.*** (2019). Comparing Arabidopsis receptor kinase and receptor protein-mediated immune signaling reveals BIK1-dependent differences. *Nature Plants*, 5(2):176–187. [DOI: 10.1038/s41477-018-0339-8](https://doi.org/10.1038/s41477-018-0339-8).
2. **Love MI, Huber W, Anders S.** (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12):550. [DOI: 10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
3. **Zhu A, Ibrahim JG, Love MI.** (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*, 35(12):2084–2092. [DOI: 10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895).
4. **Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, *et al.*** (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2(3):100141. [DOI: 10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141).
5. **Szklarczyk D, Kirsch R, Koutrouli M, Nastou K, Mehryary F, Hachilif R, Gable AL, Fang T, Doncheva NT, Pyysalo S, *et al.*** (2023). The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research*, 51(D1):D638–D646. [DOI: 10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000).

---

## Lisensi & Penulis
- **Lisensi:** Repositori ini didistribusikan di bawah lisensi [MIT License](LICENSE).
- **Penulis:** Hanif Saputra (2026).
