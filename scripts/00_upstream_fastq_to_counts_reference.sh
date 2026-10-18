#!/usr/bin/env bash
# ==============================================================================
# PIPELINE HULU RNA-SEQ (TAHAP 1 s.d. TAHAP 4) - SKRIP REFERENSI REPRODUKSIBEL
# Dataset   : GSE133053 / BioProject PRJNA549850 / SRA SRP201971 (Wan et al., 2019)
# Organisme : Arabidopsis thaliana (Col-0, TAIR10 / Araport11)
# Kontras   : Perlakuan Elisitor Kitin 1 jam (C6_1h) vs Kontrol Air 1 jam (H2O_1h)
# Catatan   : Pada lingkungan Windows lokal (eksekusi dimulai dari matriks count
#             resmi NCBI GEO GSE133053_RAW.tar), skrip bash ini disediakan untuk
#             reproduksibilitas penuh Tahap 1-4 dari FASTQ mentah di Linux/WSL.
# ==============================================================================

set -euo pipefail

# Konfigurasi resource komputasi (Batas komputasi: 8 GB RAM, 4 core)
THREADS=4
SEED=42

mkdir -p data/fastq_raw data/fastq_trimmed data/reference data/bam \
         qc/fastqc_raw qc/multiqc_raw qc/fastqc_trimmed qc/multiqc_trimmed \
         results logs

# Daftar 8 sampel fokus (1 jam Kitin C6 vs 1 jam Kontrol H2O, masing-masing 4 replikat biologis)
# SRR diperoleh langsung dari SRA Run Selector (PRJNA549850 / metadata.csv)
SRR_LIST=(
  "SRR9331441" # GSM3899066 - H2O_1h_Rep1
  "SRR9331400" # GSM3899024 - H2O_1h_Rep2
  "SRR9331422" # GSM3899047 - H2O_1h_Rep3
  "SRR9331427" # GSM3899052 - H2O_1h_Rep4
  "SRR9331398" # GSM3899022 - C6_1h_Rep1
  "SRR9331409" # GSM3899033 - C6_1h_Rep2
  "SRR9331412" # GSM3899036 - C6_1h_Rep3
  "SRR9331431" # GSM3899056 - C6_1h_Rep4
)

# ------------------------------------------------------------------------------
# TAHAP 1: Unduh Data FASTQ dari NCBI SRA (SRA Toolkit v3.0.10)
# ------------------------------------------------------------------------------
echo "[TAHAP 1] Mengunduh data SRA dengan prefetch dan fasterq-dump..."
for SRR in "${SRR_LIST[@]}"; do
  prefetch "${SRR}" --output-directory data/sra_cache
  fasterq-dump "data/sra_cache/${SRR}/${SRR}.sra" \
    --outdir data/fastq_raw \
    --threads "${THREADS}" \
    --split-files \
    --progress
  gzip -f "data/fastq_raw/${SRR}"*.fastq
done

# ------------------------------------------------------------------------------
# TAHAP 2: Quality Control (QC) Read Mentah (FastQC v0.12.1 & MultiQC v1.21)
# ------------------------------------------------------------------------------
echo "[TAHAP 2] Menjalankan FastQC dan MultiQC pada FASTQ mentah..."
fastqc data/fastq_raw/*.fastq.gz \
  --outdir qc/fastqc_raw \
  --threads "${THREADS}"

multiqc qc/fastqc_raw \
  --outdir qc/multiqc_raw \
  --filename multiqc_raw_report.html

# ------------------------------------------------------------------------------
# TAHAP 3: Trimming Adapter dan Kualitas Rendah (fastp v0.23.4)
# Parameter:
#   -q 20          : Phred quality minimal Q20 (akurasi base call >= 99%)
#   -u 30          : Buang read jika >30% basa memiliki kualitas < Q20
#   -l 36          : Buang read yang panjangnya < 36 bp pasca-trimming untuk
#                    mencegah multi-mapping ambigu pada genom TAIR10
#   --cut_tail     : Sliding window 4 bp dari ujung 3' (mean quality < Q20)
# ------------------------------------------------------------------------------
echo "[TAHAP 3] Menjalankan trimming dengan fastp dan QC ulang..."
for SRR in "${SRR_LIST[@]}"; do
  fastp \
    -i "data/fastq_raw/${SRR}.fastq.gz" \
    -o "data/fastq_trimmed/${SRR}_trimmed.fastq.gz" \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 30 \
    --length_required 36 \
    --cut_tail \
    --cut_tail_window_size 4 \
    --cut_tail_mean_quality 20 \
    --thread "${THREADS}" \
    --html "qc/fastqc_trimmed/${SRR}_fastp.html" \
    --json "qc/fastqc_trimmed/${SRR}_fastp.json" \
    2> "logs/${SRR}_fastp.log"
done

fastqc data/fastq_trimmed/*_trimmed.fastq.gz \
  --outdir qc/fastqc_trimmed \
  --threads "${THREADS}"

multiqc qc/fastqc_trimmed \
  --outdir qc/multiqc_trimmed \
  --filename multiqc_trimmed_report.html

# ------------------------------------------------------------------------------
# TAHAP 4: Alignment dan Kuantifikasi Gen
# Referensi Resmi: Ensembl Plants r58 / TAIR10 (Arabidopsis thaliana Col-0)
# Jalur A (Standar Genomik): HISAT2 v2.2.1 -> samtools v1.19 -> featureCounts (Subread v2.0.6)
# Catatan: Dalam publikasi asli GSE133053 (Wan et al., 2019), penulis memetakan
# read ke transkrip representatif TAIR10 Col-0 menggunakan Bowtie2 + eXpress
# dan menghasilkan matriks unique counts (GSE133053_RAW.tar).
# ------------------------------------------------------------------------------
echo "[TAHAP 4] Mengunduh genom referensi TAIR10 dan anotasi GTF dari Ensembl Plants..."
wget -c -O data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
wget -c -O data/reference/Arabidopsis_thaliana.TAIR10.58.gtf.gz \
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz"

gunzip -fk data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip -fk data/reference/Arabidopsis_thaliana.TAIR10.58.gtf.gz

hisat2-build -p "${THREADS}" \
  data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  data/reference/TAIR10_index

for SRR in "${SRR_LIST[@]}"; do
  hisat2 -p "${THREADS}" --seed "${SEED}" --dta \
    -x data/reference/TAIR10_index \
    -U "data/fastq_trimmed/${SRR}_trimmed.fastq.gz" \
    --summary-file "logs/${SRR}_hisat2_summary.txt" \
    | samtools sort -@ "${THREADS}" -o "data/bam/${SRR}.sorted.bam" -
  samtools index "data/bam/${SRR}.sorted.bam"
done

# Kuantifikasi tingkat gen dengan featureCounts (Illumina TruSeq unstranded: -s 0)
featureCounts -T "${THREADS}" \
  -s 0 \
  -t exon \
  -g gene_id \
  -a data/reference/Arabidopsis_thaliana.TAIR10.58.gtf \
  -o results/featureCounts_raw_matrix.txt \
  data/bam/*.sorted.bam \
  2> logs/featureCounts.log

echo "[SELESAI] Matriks count tersimpan di results/featureCounts_raw_matrix.txt"
