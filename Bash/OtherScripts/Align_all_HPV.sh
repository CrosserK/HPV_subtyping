#!/bin/bash
set -e
set -u
#set -x
set -o pipefail

#Choose FASTQ file(s)
FASTQList="Pt_48.IonXpress_054 Pt_97.IonXpress_062 Pt_1028.IonXpress_082"

#Choose References & if they should be indexed
All_references=FALSE # Use all references in reference folder
SpecificReferences="K02718.1" # Only used if above line is All_references=FALSE. Write 1 if already aligned or name
Index_references=FALSE

# Include alignment or define specific bam files
Align=TRUE # If false, will only do samtools depth and flagstats. Below params will have to be defined. 
SpfcBamF=Fra_Cecilie #Folder. Only used if Align=FALSE.
SpfcBamFiles="Pt_10.IonXpress_038 Pt_11_DNA.IonXpress_085" #Files. Only used if Align=FALSE.
Name_file_with_reference=FALSE # Define if already aligned files should have their references named in filename. Only used if Align=FALSE.
Reference_name="hg19_HPV_all" # Only used if Align=FALSE.

# Include variant calling


#Define main folder
TargetFolder='/home/pato/Skrivebord/HPV16_projekt'

# Generate folder structure
# mkdir -p $TargetFolder/{Aligned/{Samfiles,BamFiles/unsorted},References,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics}} 
# -p option makes sure that it will not delete or change existing folders

cd $TargetFolder

#Define Folders
AlaF=Aligned;RefF=References;RefdF=ReferenceDetails;SeqF=FASTQ;AnaF=Analysis;QualF=$AnaF/Qual;DepthF=$AnaF/Depth
FlagF=$AnaF/Flagstats;DupF=$AnaF/DuplicateMetrics;SamF=$AlaF/Samfiles;BamF=$AlaF/BamFiles

# Get names of all references
ls $RefF > $RefdF/HPV_subtypeList.txt
awk 'NR % 7 == 1' $RefdF/HPV_subtypeList.txt > $RefdF/HPV_subtypeList_temp1.txt
sed 's/.fasta//g' < $RefdF/HPV_subtypeList_temp1.txt > $RefdF/HPV_subtypeList.txt
rm $RefdF/HPV_subtypeList_temp1.txt

#mv $HPV_subtypeList_temp1.txt $RefdF/HPV_subtypeList.txt

# Define if using all references in folder or only specific ones
if [ $All_references == TRUE ]
then
	RefList=$(< $RefdF/HPV_subtypeList.txt)

else
	RefList=$SpecificReferences
fi

#Indexing references
if [ $Index_references == TRUE ]
then
	for refType in $RefList; do
	bwa index -a bwtsw $RefF/${refType}.fasta
	samtools faidx $RefF/${refType}.fasta
	done
fi

DepthFile=Depths_$(date +%H-%M-%S_%d-%m-%Y)
touch $DepthF/$DepthFile

#Testparams:###########################
# FASTQList="Pt_48.IonXpress_054 Pt_97.IonXpress_062 Pt_1028.IonXpress_082"
# TargetFolder='/home/pato/Skrivebord/HPV16_projekt'
# cd $TargetFolder
# Align=TRUE
# FASTQList="Pt_48.IonXpress_054"
# SpecificReferences="AF402678.1"
#Husk at enable bwa mem igen
# Enable samtools igen
#######################################

# Running main program
if [ $Align == TRUE ]
then
for FASTQ in $FASTQList; do
    for refType in $RefList; do 
    stdBamOut=$BamF/${FASTQ}_${refType}_Sorted.bam
    bwa mem $RefF/${refType}.fasta $SeqF/${FASTQ}.fastq > $SamF/${FASTQ}_$refType.sam ;
    samtools view -Sb $SamF/${FASTQ}_$refType.sam > $BamF/unsorted/${FASTQ}_$refType.bam ; 
    samtools sort $BamF/unsorted/${FASTQ}_$refType.bam -o $stdBamOut
    #Qual data
    samtools flagstat $BamF/${FASTQ}_${refType}_Sorted.bam -@ 12 > $FlagF/${FASTQ}_${refType}.txt
    DepVal=$(samtools depth $BamF/${FASTQ}_${refType}_Sorted.bam | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
    echo -e $FASTQ'\t'$refType'\t'"Depth:"'\t'$DepVal >> $DepthF/$DepthFile  

    #Duplicate markup
    java -Xmx28G -jar ~/picard.jar MarkDuplicates I=$stdBamOut \
    O=$BamF/${FASTQ}_${refType}_Sorted_DMed.bam M=$DupF/${FASTQ}_${refType}_DupMetrics.txt
    rm $stdBamOut

    stdBamOut=$BamF/${FASTQ}_${refType}_Sorted_DMed.bam

    samtools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --output-tags DP,AD -f $RefF/${refType}.fasta --BCF $stdBamOut | \
    bcftools call --consensus-caller --variants-only --pval-threshold 1.0 -Ob > $BamF/${FASTQ}_${refType}DMed.bcf

    done
done
else

    if [[ $Name_file_with_reference == TRUE ]]; then
        RefList=$Reference_name
    else
        RefList="1"
    fi

# Script without alignment on specifically named files
# Set bam subfolder
BamF=$AlaF/BamFiles/$SpfcBamF

for BAM in $SpfcBamFiles; do
    for refType in $RefList; do
        samtools flagstat $BamF/${BAM}.bam -@ 12 > $FlagF/${BAM}_${refType}.txt
        DepVal=$(samtools depth $BamF/${BAM}.bam | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }')
        echo -e $BAM'\t '$refType'\t '"Depth:"'\t'$DepVal >> $DepthF/$DepthFile

    done
done
fi


# Get number of mismatches
cat '/home/pato/Skrivebord/HPV16_projekt/VariantCalls/TSVC_variants_IonXpress_058.vcf' | grep -v "#" | grep "" -c



