#!/bin/bash
set -e
set -u # Exit if undeclared variables

# define help function
function help ()
{
	printf "Usage: %s: some var is missing 
	[-h|--help]\n" "$(basename "$0")" >&2
    exit
}

# Reading input parameters
while [[ "$#" -gt 0 ]]; do
	case "$1" in
		-f|--fastq) FastQFile="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-s|--TopRunName) TopRunName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -m|--mainf) MainF="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -r|--reference) Reference="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-t|--typetier) typeTier="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -b|--bedfile) BedFolder="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -n|--bedpoolN) BedPoolN="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

FQName=$(basename ${FastQFile} .fastq)


# Define locations
QCF=$MainF/QC; 
RefF=$MainF/References; 
LogF=$QCF/Logs;
IdxF=$QCF/Indexing
FlagF=$QCF/Flagstats;
ResultsF=$MainF/Results/$TopRunName	

###################### Find genotypecall ref ##########################

echo -e "Using reference" "$Reference" for $FQName
echo [$(date +"%d-%m-%Y %H:%M:%S")] "Using reference" "$Reference" for $FQName >> $LogF/${TopRunName}.txt

mkdir -p $ResultsF/$FQName
workD=$ResultsF/$FQName
mkdir -p $workD/$Reference
currentF=$workD/$Reference
mkdir -p "$currentF"/ResultFiles

# Get number of amplicon pools 
readarray -t amplPools < <(find $BedFolder -maxdepth 1 -name "pool*.bed")
BedPoolN=${#amplPools[@]}

if [ $BedPoolN -gt 1 ];
then
    # Using bedfile to cut primer sequences
    for i in 1 $BedPoolN; 
    do

        echo Splitting in "$BedPoolN" pools

        # Choose file 
        tmpBedFile=$FQF/"pool"${i}".bed"

        # Splitting to pools. This is done in order to allow better clipping where reads going into a neighbour primer are kept. 
        # -F is minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
        poolWithOneRefOnly=${tmpBedFile%.bed}_${Reference}.bed
        grep "$Reference" $tmpBedFile > $poolWithOneRefOnly

        bedtools intersect -a $BamFile -b $poolWithOneRefOnly > "${BamFile%.bam}"_intersect${i}.bam -F 1.0
        BamInt="${BamFile%.bam}"_intersect${i}.bam

        clipped_out=${BamInt%.bam}_clip.bam
        # Now clip with samtools ampliconclip
        samtools ampliconclip -o $clipped_out -b $tmpBedFile $BamInt --hard-clip --both-ends --filter-len 80 --tolerance 0
        # Sorting by name for samtools fixmate
        samtools sort -n $clipped_out > ${clipped_out%.bam}.sort.bam 
        # Fixing new TLEN tags
        samtools fixmate ${clipped_out%.bam}.sort.bam $clipped_out
        # Fixing new MD and NM aux tags
        samtools calmd $clipped_out $Ref_FASTA > ${clipped_out%.bam}_calmd.bam
        # BGZIP for samtools sort
        bgzip ${clipped_out%.bam}_calmd.bam
        # Sorting by coordinates
        samtools sort ${clipped_out%.bam}_calmd.bam.gz > ${clipped_out%.bam}_calmd.sort.bam.gz
        # Indexing
        samtools index ${clipped_out%.bam}_calmd.sort.bam.gz

        clipped_out=${clipped_out%.bam}_calmd.sort.bam.gz

        # HaplotypeCaller nødvendigheder
        # Tilføjer tags, nødvendig for GATK, da der ikke arbejdes med uBAM filer 
        # RGLB = Read group library, RGPL = platform, RGPU = platform unit, RGSM = sample name
        java -jar picard.jar AddOrReplaceReadGroups \
        I=$clipped_out \
        O="${clipped_out%.bam}".readGroupFix.bam \
        RGLB=lib1 \
        RGPL=IonTorrent \
        RGPU=unit1 \
        RGSM=1

        BamOut="${clipped_out%.bam}".readGroupFix.bam

        samtools index $BamOut

        done

        # Does not merge bam files from each pool as cuplicates can occur, because
        # both splits can contain some of the same reads
else

    Ref_FASTA=$RefF/IndexedRef/"${Reference}"/"${Reference}".fasta # Find reference for picard
    echo Aligning with reference ${Ref_FASTA##*/}
    echo [$(date +"%d-%m-%Y %H:%M:%S")] Aligning with reference "${Ref_FASTA##*/}" >> $LogF/${TopRunName}.txt
    BamFile=$currentF/${Reference}.bam

    # Aligner
    bwa mem -t 12 -v 2 -M "$Ref_FASTA" $FastQFile > "${BamFile%bam}"sam # -v = verbosity, 2 for errors og warnings kun
    # Sort bamfile
    samtools sort "${BamFile%bam}"sam -o "${BamFile%bam}"sort.bam 

    java -Xmx28G -jar ~/picard.jar MarkDuplicates -I "${BamFile%bam}"sort.bam -O "${BamFile%bam}"sort.dup.bam \
    -M "${BamFile%.bam}"_DupMetrics.txt -REMOVE_DUPLICATES false 2> $IdxF/markdup_log_"${Reference}".txt

    java -Xmx28G -jar ~/picard.jar MarkDuplicates -I "${BamFile%bam}"sort.bam -O "${BamFile%bam}"sort.dup_rm.bam \
    -M "${BamFile%.bam}"_DupMetrics.txt -REMOVE_DUPLICATES true 2> $IdxF/markdup_log_"${Reference}".txt

    # Rename & index nodupmarked file
    samtools index "${BamFile%bam}"sort.bam 
    samtools index "${BamFile%bam}"sort.dup.bam
    samtools index "${BamFile%bam}"sort.dup_rm.bam

    samtools flagstats "${BamFile%bam}"sort.bam -O tsv > $FlagF/${FQName}.tsv

    # Her kan angives om dup eller ikke dup skal bruges til resten af kørsel (.dup.)
    BamFile="${BamFile%bam}"sort.bam

    cp $BamFile $currentF/ResultFiles/${FQName}_${Reference}.sort.bam;
    cp ${BamFile}.bai $currentF/ResultFiles/${FQName}_${Reference}.sort.bam.bai;

fi