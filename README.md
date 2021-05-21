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
8. (In progress...) Call variants with GATK HaplotypeCaller
9. (In progress...) Filter variants with GATK VariantRecalibrator in SNP mode + ApplyRecalibration
10. (In progress...) Filter variants again but this time with GATK VariantRecalibrator in indel mode + ApplyRecalibration
11. (In progress...) Remove filtered variants from vcf with GATK SelectVariants --excludeFiltered
12. (In progress...), optional) Evaluate variantcalling against similar data with GATK VariantEval
13. Find number of mismatches in vcf files
14. Find most likely subtype with sigma (sigma-build + sigma solve) & compare with number of mismatches from vcf
15. Find amino acid changes by variants in the subtypes in R
    - Prepare gff3 file with GenomicFeatures library -> makeTxDbFromGFF
    - VariantAnnotation library -> predictCoding
16. (In progress...) Create html report of sequenced subtype'
