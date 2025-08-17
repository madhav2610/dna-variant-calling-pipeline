#!/bin/bash

# ===== STEP 10: Sort BAM file by name for FixMate =====
samtools sort -n -o aligned_name_sorted.bam aligned.bam

# ===== STEP 11: Fix mate information =====
samtools fixmate aligned_name_sorted.bam fixed.bam

# ===== STEP 12: Sort BAM file by coordinate =====
samtools sort -o aligned_coord_sorted.bam fixed.bam

# ===== STEP 13: Mark duplicates using GATK =====
gatk MarkDuplicates \
  -I aligned_coord_sorted.bam \
  -O dedup.bam \
  -M dedup_metrics.txt

# ===== STEP 14: Index BAM file =====
samtools index dedup.bam

# ===== STEP 15: Create sequence dictionary for reference genome =====
gatk CreateSequenceDictionary -R chr20.fasta

# ===== STEP 16: Index reference genome =====
samtools faidx chr20.fasta
