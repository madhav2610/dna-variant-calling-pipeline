#!/bin/bash

# ===== STEP 21: Annotate variants (add rsIDs) =====
gatk VariantAnnotator \
  -R chr20.fasta \
  -V raw_variants.vcf \
  -O annotated_variants.vcf \
  --dbsnp dbsnp_chr20.vcf.gz

# ===== STEP 22: Separate SNPs and INDELs =====
gatk SelectVariants \
  -R chr20.fasta \
  -V annotated_variants.vcf \
  --select-type-to-include SNP \
  -O raw_snps.vcf

gatk SelectVariants \
  -R chr20.fasta \
  -V annotated_variants.vcf \
  --select-type-to-include INDEL \
  -O raw_indels.vcf

# ===== STEP 23: Filter SNPs =====
gatk VariantFiltration \
  -R chr20.fasta \
  -V raw_snps.vcf \
  -O filtered_snps.vcf \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "snp_filter"

# ===== STEP 24: Filter INDELs =====
gatk VariantFiltration \
  -R chr20.fasta \
  -V raw_indels.vcf \
  -O filtered_indels.vcf \
  --filter-expression "QD < 2.0 || FS > 200.0" \
  --filter-name "indel_filter"

# ===== STEP 25: Count SNPs and INDELs before and after filtering =====
grep -v "^#" raw_snps.vcf | wc -l
grep -v "^#" filtered_snps.vcf | wc -l

grep -v "^#" raw_indels.vcf | wc -l
grep -v "^#" filtered_indels.vcf | wc -l
