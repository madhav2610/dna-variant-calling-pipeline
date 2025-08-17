#!/bin/bash

# ===== STEP 4: Run FastQC on FASTQ files =====
fastqc SRR12023506_1.fastq
fastqc SRR12023506_2.fastq

# ===== STEP 5: Trim adapters and low-quality bases using Trimmomatic =====
trimmomatic PE -threads 8 \
  SRR12023506_1.fastq SRR12023506_2.fastq \
  SRR12023506_1_paired.fastq SRR12023506_1_unpaired.fastq \
  SRR12023506_2_paired.fastq SRR12023506_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# ===== STEP 6: Repair read pairs using BBMap repair.sh =====
repair.sh \
  in1=SRR12023506_1_paired.fastq \
  in2=SRR12023506_2_paired.fastq \
  out1=SRR12023506_1_repaired.fastq \
  out2=SRR12023506_2_repaired.fastq \
  outs=singles.fastq

