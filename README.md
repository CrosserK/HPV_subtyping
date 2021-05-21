# HPV_subtyping (In development)
Workflow to get from FASTQ files to calling the HPV16 subtype by

### Workflow

1. Trimming low qual bases from ends and reads <75 and >150 bp with cutadapt
2. (optional) Quality checking with FastQC
3. Indexing all given reference subtypes with sigma-index
4. Aligning all fastq files with sigma-align (uses Bowtie2 for alignment)
5. Sort bamfiles with samtools sort
6. Mark duplicates with picard tools
7. Index bam with samtools index
8. (not done yet) Call variants with GATK HaplotypeCaller
9. (not done yet, only cohort data?) GATK PairHMM
10. (not done yet) Filter variants with GATK VariantRecalibrator in SNP mode + ApplyRecalibration
11. (not done yet) Filter variants again but this time with GATK VariantRecalibrator in indel mode + ApplyRecalibration
12. Remove filtered variants from vcf with GATK SelectVariants --excludeFiltered
13. (not done yet, optional) Evaluate variantcalling against similar data with GATK VariantEval
14. Find number of mismatches in vcf files
15. Find most likely subtype with sigma (sigma-build + sigma solve) & compare with number of mismatches from vcf
16. Find amino acid changes by variants in the subtypes in R
    - Prepare gff3 file with GenomicFeatures library -> makeTxDbFromGFF
    - VariantAnnotation library -> predictCoding
17. (not done yet) Create html report of sequenced subtype
