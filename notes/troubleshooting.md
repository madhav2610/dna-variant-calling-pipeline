# Troubleshooting Guide

This file documents the common errors I faced while building this pipeline and their fixes.  
If you encounter the same issues, check here first.

---

## 1. BWA Alignment Error (reads differently paired)

**Error:**  
BWA stopped after ~17 GB with error about reads being differently paired.

**Cause:**  
FastQC showed adapter content was near zero, but the read counts in R1 and R2 files were different.

**Fix:**  
Used `repair.sh` from BBMap to fix mismatched pairs:  
```bash
repair.sh in1=SRR12023506_1.fastq in2=SRR12023506_2.fastq out1=SRR12023506_1_repaired.fastq out2=SRR12023506_2_repaired.fastq outs=singles.fastq
```

 ## 2.  BQSR Error: Contig Name Mismatch

**Error:**
A USER ERROR has occurred: Input files reference and known sites have incompatible contigs


**Cause:**
Reference contigs had names like CM000682.2 while dbSNP used 1, 2, 3 ... X, Y.

**Fix:**  
Rename contigs in dbSNP file using bcftools:
```bash
bcftools annotate --rename-chrs rename.txt input.vcf.gz -Oz -o fixed.vcf.gz
tabix -p vcf fixed.vcf.gz
(where rename.txt maps contig IDs).
```
 ## 3. GATK Error: "Read group platform (PL) is missing"

**Error:**
A USER ERROR has occurred: Read group platform (PL) not set in BAM header


**Cause:**
BAM file did not contain @RG header with platform = ILLUMINA.

**Fix:**
Add read groups using Picard:
```bash
picard AddOrReplaceReadGroups \
  I=aligned.bam \
  O=rg_added.bam \
  RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=sample1
```
 ## 4. GATK Installation / Path Issues

**Error:**
command not found: gatk or older version running.

**Fix:**
Installed miniconda

Created environment:
```bash
conda create -n bioinfo
conda activate bioinfo
conda install -c bioconda gatk4 picard samtools bcftools fastqc trimmomatic bbmap
```
Added environment to PATH.

## 5. Memory Issues with SRA Toolkit

**Error:**
fasterq-dump crashing on large SRA files.

**Fix:**

Use --split-files and -e <threads> to speed up.

Ensure enough disk space (>100 GB free).