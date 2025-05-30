#!/bin/bash

# Kallisto RNA-Seq Pipeline (Clean Version)
# Author: Menatallah Adel Abdelmagid
# System used: Ubuntu 20.04 (via VirtualBox) on HP laptop


# Step 1: Create project directory and navigate into it
mkdir transcriptome 

cd transcriptome

# Step 2: Download and extract human transcriptome reference (Ensembl GRCh38)
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

# Step 3: Install kallisto (using conda)
conda install -c bioconda kallisto

# Step 4: Build kallisto index
kallisto index -i human2.idx Homo_sapiens.GRCh38.cdna.all.fa

# Step 5: Quality control on raw FASTQ files
fastqc *.fastq

# Step 6: Trim adapters and low-quality bases using Trimmomatic
trimmomatic PE SRR1000_R1.fastq SRR1000_R2.fastq \
SRR1000_R1_trimmed.fq.gz SRR1000_R1_drop.fq.gz \
SRR1000_R2_trimmed.fq.gz SRR1000_R2_drop.fq.gz \
ILLUMINACLIP:./adaptor.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100

# Step 7: Quality control on trimmed files
mkdir fastqc_trimmedreports
for file in *_trimmed.fq.gz; do
    fastqc "$file" -o ./fastqc_trimmedreports/
done

# Step 8: Move index and reference file to working folder (optional)
mv human2.idx Homo_sapiens.GRCh38.cdna.all.fa ./Transcriptome/

# Step 9: Run kallisto quantification
kallisto quant -i human2.idx -o SRR1000 SRR1000_R1_trimmed.fq.gz SRR1000_R2_trimmed.fq.gz
