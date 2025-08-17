#!/bin/bash

# ===== STEP 7: Index the reference genome with BWA =====
bwa index chr20.fasta

# ===== STEP 8: Align reads to reference using BWA =====
bwa mem -t 8 chr20.fasta \
  SRR12023506_1_repaired.fastq \
  SRR12023506_2_repaired.fastq > aligned.sam

# ===== STEP 9: Convert SAM to BAM =====
samtools view -bS aligned.sam > aligned.bam
