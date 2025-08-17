#!/bin/bash

# ===== STEP 17: Download dbSNP VCF =====
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GCF_000001405.39.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GCF_000001405.39.gz.tbi

# ===== STEP 18: Extract chr20 from dbSNP file =====
bcftools view -r chr20 GCF_000001405.39.gz -Oz -o dbsnp_chr20.vcf.gz
tabix -p vcf dbsnp_chr20.vcf.gz

# ===== STEP 19: Run Base Quality Score Recalibration (BQSR) =====
gatk BaseRecalibrator \
  -R chr20.fasta \
  -I dedup.bam \
  --known-sites dbsnp_chr20.vcf.gz \
  -O recal_data.table

gatk ApplyBQSR \
  -R chr20.fasta \
  -I dedup.bam \
  --bqsr-recal-file recal_data.table \
  -O recal.bam

# ===== STEP 20: Call variants using HaplotypeCaller =====
gatk HaplotypeCaller \
  -R chr20.fasta \
  -I recal.bam \
  -O raw_variants.vcf
