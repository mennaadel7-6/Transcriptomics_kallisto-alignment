# Transcriptomics_kallisto-alignment

This repository contains a streamlined RNA-seq data processing pipeline using Kallisto for transcript quantification. The pipeline covers downloading data, quality control, trimming, and quantification steps using standard bioinformatics tools on a Linux system.

---

## Overview

This pipeline performs the following main steps:

- Download and prepare reference transcriptome  
- Quality control of raw FASTQ files using FastQC and MultiQC  
- Trim adapters and low-quality sequences with Trimmomatic  
- Quality control on trimmed FASTQ files  
- Quantify transcript abundance using Kallisto  
- Merge individual sample quantifications into a single count matrix  

---

## Main Commands Summary

```bash
# 1. Download and Prepare Reference
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
kallisto index -i human2.idx Homo_sapiens.GRCh38.cdna.all.fa

# 2. Download Raw RNA-seq Data and Create FASTQ files
cat srr_list.txt | xargs -I {} fasterq-dump --split-files {}
gzip *.fastq

# 3. Quality Control on Raw Data
fastqc *.fastq.gz
multiqc .

# 4. Trim Low-Quality Reads and Adapters (Trimmomatic)
for R1 in *_1.fastq.gz; do
  SAMPLE="${R1%_1.fastq.gz}"
  R2="${SAMPLE}_2.fastq.gz"
  trimmomatic PE "$R1" "$R2" \
    "${SAMPLE}_R1_trimmed.fq.gz" "${SAMPLE}_R1_drop.fq.gz" \
    "${SAMPLE}_R2_trimmed.fq.gz" "${SAMPLE}_R2_drop.fq.gz" \
    ILLUMINACLIP:./adaptor.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100
done

# 5. Quality Control on Trimmed Data
fastqc *_trimmed.fq.gz
multiqc .

# 6. Run Kallisto Quantification
for R1 in *_R1_trimmed.fq.gz; do
  R2="${R1%_R1_trimmed.fq.gz}_R2_trimmed.fq.gz"
  sample="${R1%_R1_trimmed.fq.gz}"
  kallisto quant -i human2.idx -o "$sample" "$R1" "$R2"
done

# 7. Merge Quantification Results into a Count Matrix
echo -ne "target_id" > count_matrix.tsv
for i in SRR*; do echo -ne "\t$i"; done >> count_matrix.tsv
echo "" >> count_matrix.tsv

awk -F'\t' '{print $1}' $(ls SRR*/abundance.tsv | head -n 1) > temp_ids.tsv
mv temp_ids.tsv count_matrix.tsv

for i in SRR*; do
  echo "$i" > temp_counts.tsv
  awk -F'\t' 'NR > 1 {print $4}' "$i/abundance.tsv" >> temp_counts.tsv
  paste count_matrix.tsv temp_counts.tsv > temp_matrix.tsv
  mv temp_matrix.tsv count_matrix.tsv
done

rm -f temp_counts.tsv temp_matrix.tsv

