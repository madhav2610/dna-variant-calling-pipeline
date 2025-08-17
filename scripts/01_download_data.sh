#!/bin/bash

# ===== STEP 1: Download reference genome (chr20) from NCBI =====
wget https://www.ncbi.nlm.nih.gov/nucleotide/CM000682.2?report=fasta&log$=seqview -O chr20.fasta

# ===== STEP 2: Download SRA sample from NCBI =====
# First install SRA Toolkit (if not installed already)
# Download from GitHub, extract, and set path
# (these steps you already did manually, so I’ll only keep main commands here)

# Download sample SRA file using prefetch
prefetch SRR12023506

# ===== STEP 3: Convert .sra to FASTQ =====
fasterq-dump SRR12023506 --split-files

# After this you should have:
# SRR12023506_1.fastq
# SRR12023506_2.fastq
