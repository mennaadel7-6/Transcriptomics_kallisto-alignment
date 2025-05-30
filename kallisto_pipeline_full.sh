#!/bin/bash

# Full Kallisto RNA-seq Workflow
# Author: Menatallah Adel Abdelmagid (this script was a group work with NGS group)
# System used: Ubuntu 20.04  on HP laptop

# -----------------------------
# 1. Download and convert SRA to FASTQ
# -----------------------------
cat srr_list.txt | xargs -I {} fasterq-dump --split-files {}

# -----------------------------
# 2. Compress FASTQ files
# -----------------------------
gzip *.fastq

# -----------------------------
# 3. Quality check using FastQC
# -----------------------------
fastqc *.fastq.gz

# -----------------------------
# 4. Generate MultiQC summary
# -----------------------------
multiqc .

# -----------------------------
# 5. Trimming using Trimmomatic
# -----------------------------
for R1 in *_1.fastq.gz; do
    SAMPLE="${R1%_1.fastq.gz}"
    R2="${SAMPLE}_2.fastq.gz"
    trimmomatic PE "$R1" "$R2" \
    "${SAMPLE}_R1_trimmed.fq.gz" "${SAMPLE}_R1_drop.fq.gz" \
    "${SAMPLE}_R2_trimmed.fq.gz" "${SAMPLE}_R2_drop.fq.gz" \
    ILLUMINACLIP:./adaptor.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100
done

# -----------------------------
# 6. FastQC After Trimming
# -----------------------------
fastqc *_trimmed.fq.gz

# -----------------------------
# 7. Run MultiQC again
# -----------------------------
multiqc .

# -----------------------------
# 8. Pseudoalignment using Kallisto
# -----------------------------
for R1 in *_R1_trimmed.fq.gz; do
    R2="${R1%_R1_trimmed.fq.gz}_R2_trimmed.fq.gz"
    sample="${R1%_R1_trimmed.fq.gz}"
    kallisto quant -i human2.idx -o "$sample" "$R1" "$R2"
done

# -----------------------------
# 9. Create Count Matrix from Kallisto outputs
# -----------------------------

# Create header row
echo -ne "target_id" > count_matrix.tsv
for i in SRR*; do echo -ne "\t$i"; done >> count_matrix.tsv
echo "" >> count_matrix.tsv

# Extract gene IDs
awk -F'\t' '{print $1}' $(ls SRR*/abundance.tsv | head -n 1) > temp_ids.tsv
mv temp_ids.tsv count_matrix.tsv

# Append TPM values from each sample
for i in SRR*; do
    echo "$i" > temp_counts.tsv
    awk -F'\t' 'NR > 1 {print $4}' "$i/abundance.tsv" >> temp_counts.tsv
    paste count_matrix.tsv temp_counts.tsv > temp_matrix.tsv
    mv temp_matrix.tsv count_matrix.tsv
done

# Cleanup
rm -f temp_counts.tsv temp_matrix.tsv
