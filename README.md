# HPV_subtyping (In development)
Workflow to get from FASTQ files to calling the HPV16 subtype, given references of possible subtypes. The worfklow will also generate a report with nucleotide and amino acid changes. 

### Workflow

1. Trimming low qual bases from ends and reads \<x and \>y bp with cutadapt
2. Quality checking with FastQC
3. Indexing all given reference subtypes with sigma-index
4. Aligning all fastq files with sigma-align (uses Bowtie2 for alignment)
5. Sort bamfiles with samtools sort
6. Mark duplicates with picard tools
7. Index bam with samtools index
8. Call variants with GATK HaplotypeCaller
9. Get vcf stats DP, QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum with sed for R graphing (NOTE: the latter 2 are not grabbed correctly with sed)
10. Hard filter variants with GATK VariantFiltration according to chosen filters (optionally use graphs)
12. Remove filtered variants from vcf with GATK SelectVariants --exclude-filtered
14. Find number of mismatches in vcf files
15. Find most likely subtype with sigma (sigma-build + sigma solve) & compare with number of mismatches from vcf
16. Find amino acid changes by variants in the subtypes in R
    - Prepare gff3 file with GenomicFeatures library -> makeTxDbFromGFF
    - VariantAnnotation library -> predictCoding
17. (In progress...) Format output to e.g. p.V600E, with codon number and c.1983A>T with nucleotide number
18. (In progress...) Create html report of sequenced subtype
