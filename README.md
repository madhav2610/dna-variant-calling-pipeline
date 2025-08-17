# DNA Variant Calling Pipeline (Chromosome 20)

This repository contains a **bioinformatics pipeline** for calling SNPs and indels from short-read sequencing data.  
It was built for **learning purposes** and demonstrates the complete workflow from raw data to filtered variants.

---

## 🚀 Workflow Overview

**Steps:**
1. **Download reference genome** (chr20 only) from NCBI.  
   - [CM000682.2](https://www.ncbi.nlm.nih.gov/nucleotide/CM000682.2)

2. **Download sequencing data** from SRA.  
   - [SRR12023506](https://www.ncbi.nlm.nih.gov/sra/?term=SRR12023506)  
   - Convert SRA → FASTQ (via `fasterq-dump`).

3. **Quality Control & Trimming**  
   - QC with `FastQC`.  
   - Trimming with `Trimmomatic`.  
   - Checked adapters → found minimal adapter contamination.  

4. **Repair & Alignment**  
   - Used `bbmap/repair.sh` to fix mismatched read pairs.  
   - Indexed reference with `bwa`.  
   - Aligned reads (`bwa mem`) → SAM → BAM.

5. **BAM Processing**  
   - Fix mate info with `samtools fixmate`.  
   - Coordinate sorting + mark duplicates (`GATK MarkDuplicates`).  
   - Create dictionary + index for reference genome.  

6. **Base Quality Score Recalibration (BQSR)**  
   - Downloaded dbSNP database (GRCh38, b151).
   - [dbSNP](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/)  
   - Extracted chr20 variants only.  
   - Fixed contig name mismatches (ref vs dbsnp).  
   - Fixed missing read group info (`picard AddOrReplaceReadGroups`).  
   - Ran BQSR.

7. **Variant Calling**  
   - `GATK HaplotypeCaller` → raw VCF.  

8. **Annotation**  
   - Added rsIDs using `VariantAnnotator`.  

9. **Variant Selection & Filtering**  
   - Split SNPs and INDELs (`SelectVariants`).  
   - Filtered using `VariantFiltration`.  
   - **Results:**  
     - SNPs: 104,446 → 100,172 (after filtering).  
     - INDELs: 20,767 → 19,678 (after filtering).  

---

## 🛠 Tools & Dependencies

- **GATK4**
- **BWA**
- **Samtools**
- **Picard**
- **FastQC**
- **Trimmomatic**
- **bbmap**
- **bcftools**
- **SRA Toolkit**

All tools can be installed with conda:

```bash
conda env create -f environment.yml
conda activate bioinfo
````

---

## 📂 Repository Contents

* `scripts/` → step-wise pipeline scripts.
* `config/` → reference links.
* `notes/troubleshooting.md` → full error log + fixes (my learning journey).
* `results_files/` → filtered variant file for chr20.
* `workflow_diagram.png` → visual overview of pipeline.

---

## 📈 Workflow Diagram

<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/53f3569f-3e2a-4594-966e-b141c9ce28d2" />


---

## Notes
- These scripts are written based on my own workflow using `chr20` reference and sample `SRR12023506`.
- File names, reference genome, and sample IDs are **hardcoded** for demonstration purposes.
- If you want to use these scripts, please edit them according to your own file names, paths, and parameters.
- The pipeline order is general for DNA variant calling, but exact commands may vary depending on your dataset and tool versions.

---

## 📎 References

* [NCBI Reference Genome](https://www.ncbi.nlm.nih.gov/nucleotide/CM000682.2)
* [SRA Toolkit](https://github.com/ncbi/sra-tools)
* [dbSNP Build 151, GRCh38](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/)
* [GATK Documentation](https://gatk.broadinstitute.org/)

---

## 👤 Author

* Madhav Takkar
* This was a **practice project** to learn DNA variant calling workflows.
* Errors, fixes, and decisions are documented in `notes/troubleshooting.md`.

```


