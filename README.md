# HPV_subtyping (In development)
Workflow to get from FASTQ files to calling the HPV16 subtype, given references of possible subtypes. The worfklow will also generate a report with nucleotide and amino acid changes. 

### Workflow

1. Put all fasta files in References top folder (subfolders are not scanned for references)
2. If multiple fastq files are to be run, use script Multiple_fastq_files_HPV_subtyping.sh
3. Run fastqQC.sh, input fastqfile
   1. Trimming low qual bases from fastqfiles from ends and reads \<x and \>y bp with cutadapt
   2. Quality checking with FastQC -> evaluate if step 1 needs to be redone with new parameters
4. Run HPV_subtyping.sh, input fastqfilename and any runname
   1. Indexing all given reference subtypes with sigma-index
   2. Aligning all fastq files with sigma-align (uses Bowtie2 for alignment)
   3. Sort bamfiles with samtools sort
   4. Mark duplicates with picard tools
   5. Index bam with samtools index
   6. Call variants with GATK HaplotypeCaller
   7. Get vcf stats DP, QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum with sed for R graphing (NOTE: the latter 2 are not grabbed correctly with sed)
5. Evaluate variant calls with Graph_vcf_stats.R
6. Hard filter variants with vcf_filter.sh (GATK VariantFiltration) according to chosen filters
7. Remove filtered variants from vcf with vcf_filt_ex.sh (GATK SelectVariants --exclude-filtered)
8. Find number of mismatches in vcf files
9. Find most likely subtype with sigma (sigma-build + sigma solve)
10. Compare sigma results with number of mismatches from vcf
11. Find amino acid changes with annotate_vcf.R
    1. Prepare gff3 file with GenomicFeatures library -> makeTxDbFromGFF
    2. VariantAnnotation library -> predictCoding
    3. Format output to e.g. p.V600E, with codon number and c.1983A>T with nucleotide number
    4. (In progress...) Create html report of sequenced subtype
