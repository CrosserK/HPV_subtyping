HPV Subtyping (in development)

Workflow to get from FASTQ files to calling the HPV geno- and subtype, given references of possible types. The worfklow will also generate a report with nucleotide and amino acid changes.
Workflow

    Put all fasta files in References top folder (subfolders are not scanned for references)
    If multiple fastq files are to be run, use script Multiple_fastq_files_HPV_subtyping.sh
    Run fastqQC.sh, input fastqfilename (NOTE: currently HPV_subtyping.sh contains all modules from step 3 to 9)
        Trimming low qual bases from fastqfiles from ends and reads <x and >y bp with cutadapt
        Quality checking with FastQC -> evaluate if step 1 needs to be redone with new parameters
    Run HPV_subtyping.sh, input fastqfilename and any runname
        Indexing all given reference subtypes with BWA index and samtools faidx
        Aligning all fastq files with BWA-MEM
        Sort bamfiles with samtools sort
        Mark duplicates with picard tools
        Index bam with samtools index
        Call variants with GATK HaplotypeCaller
        Get vcf stats DP, QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum with sed for R graphing (NOTE: the latter 2 are not grabbed correctly with sed)
    Evaluate variant calls with Graph_vcf_stats.R
    Hard filter variants with vcf_filter.sh (GATK VariantFiltration) according to chosen filters
    Remove filtered variants from vcf with vcf_filt_ex.sh (GATK SelectVariants --exclude-filtered)
    Find number of mismatches in vcf files
    Find most likely subtype with Viral strain caller (Not chosen yet)
    Compare Viral strain caller results with number of mismatches from vcf
    Find amino acid changes with annotate_vcf.R
        Prepare gff3 file with GenomicFeatures library -> makeTxDbFromGFF
        VariantAnnotation library -> predictCoding
        Format output to e.g. p.V600E, with codon number and c.1983A>T with nucleotide number
        (In progress...) Create html report of sequenced subtype
