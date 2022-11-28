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
		-f|--fastq) FQAddr="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -m|--mainf) MainF="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-b|--bamfile) BamFile="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -r|--reference) Reference="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -a|--BedFolder) BedFolder="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

RefF=$MainF/References; 

# Get number of amplicon pools 
readarray -t amplPools < <(find $BedFolder -maxdepth 1 -name 'pool[0-9]?bed')
BedPoolN=${#amplPools[@]}

# Using bedfile to cut primer sequences
for (( i=1; i<=$BedPoolN; i++ ))
do

    if [ $BedPoolN == 2 ];
    then
        echo Splitting in "$BedPoolN" pools and removing primers
    elif [ $BedPoolN == 1 ];
    then
        echo Removing primers
    fi

    Ref_FASTA=$RefF/IndexedRef/"${Reference}"/"${Reference}".fasta # Find reference for picard

    # Choose file 
    tmpBedFile=$BedFolder/"pool"${i}".bed"

    # Splitting to pools. This is done in order to allow better clipping where reads going into a neighbour primer are kept. 
    # -F is minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).
    poolWithOneRefOnly=${tmpBedFile%.bed}_${Reference}.bed
    grep "$Reference" $tmpBedFile > $poolWithOneRefOnly

    bedtools intersect -a $BamFile -b $poolWithOneRefOnly > "${BamFile%.bam}"_intersect${i}.bam
    BamInt="${BamFile%.bam}"_intersect${i}.bam

    # Splitting to pools. This is done in order to allow better clipping where reads going into a neighbour primer are kept. 
    bedtools intersect -a $BamFile -b $tmpBedFile > "${BamFile%.bam}"_intersect${i}.bam -F 0.95
    samtools index $BamInt

    # Edit the bed file to make the positions span primers instead of amplicons
    while IFS= read -r line
    do
        chr=$(echo $line | awk '{print $1'})
        startpos=$(echo $line | awk '{print $2'})
        endpos=$(echo $line | awk '{print $3'})
        
        startfront=$((startpos-25))
        startend=$((startpos))
        endfront=$((endpos))
        endend=$((endpos+25))

        echo -e $chr"\t"$startfront"\t"$startend >> "${tmpBedFile%.bed}"_primers.bed
        echo -e $chr"\t"$endfront"\t"$endend >> "${tmpBedFile%.bed}"_primers.bed

    done < "$tmpBedFile"

    # Clean up the bam (make unmapped reads have flag 0)
    java -jar ~/picard.jar CleanSam \
    I=$BamInt \
    O=${BamInt%.bam}_clean.bam

    mv ${BamInt%.bam}_clean.bam $BamInt

    clipped_out=${BamInt%.bam}_clip.bam
    # Now clip with samtools ampliconclip
    samtools ampliconclip -o $clipped_out -b "${tmpBedFile%.bed}"_primers.bed $BamInt --hard-clip --both-ends --tolerance 0 --filter-len 50
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
    java -jar ~/picard.jar AddOrReplaceReadGroups \
    I=$clipped_out \
    O="${clipped_out%.bam}".readGroupFix.bam \
    RGLB=lib1 \
    RGPL=IonTorrent \
    RGPU=unit1 \
    RGSM=1 \
    VALIDATION_STRINGENCY=LENIENT

    BamOut="${clipped_out%.bam}".readGroupFix.bam

    samtools index $BamOut


done

if [ $BedPoolN == 2 ];
then
    samtools merge -o "${BamFile%.bam}"_clipped_merge.sort.bam.gz "${BamFile%.bam}"_intersect1_clip_calmd.sort.bam.gz "${BamFile%.bam}"_intersect2_clip_calmd.sort.bam.gz
    # Then convert back to fastq for use of the subtyping
    samtools fastq "${BamFile%.bam}"_clipped_merge.sort.bam.gz > $FQAddr
    echo Completed fastq conversion
else
    samtools fastq "${BamOut}" > $FQAddr
    echo Completed fastq conversion
fi


