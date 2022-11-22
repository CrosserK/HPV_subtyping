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
