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
		-a|--ampliconref) AmpliconRef="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -b|--bamfile) BamFile="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-t|--typetier) typeTier="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

FQName=$(basename ${FastQFile} .fastq)

# Define locations
QCF=$MainF/QC; 
LogF=$QCF/Logs
RefF=$MainF/References; 
ResultsF=$MainF/Results/$TopRunName	

################## VARIANT DISCOVERY #################
touch $ResultsF/$FQName/MismatchCounts_"${FQName}".txt
touch $ResultsF/$FQName/MismatchCounts_filt_"${FQName}".txt
MMcountFile=$ResultsF/$FQName/MismatchCounts_"${FQName}".txt
MMcountFile_filt=$ResultsF/$FQName/MismatchCounts_filt_"${FQName}".txt # Filtered MM count file

workD=$ResultsF/$FQName
currentF=$workD/"$Reference"
mkdir -p $currentF/ResultFiles/
Ref_FASTA=$RefF/IndexedRef/"${Reference}"/"${Reference}".fasta # Find reference for picard
echo Calling variants with reference ${Ref_FASTA##*/}
echo [$(date +"%d-%m-%Y %H:%M:%S")] Calling variants with reference ${Ref_FASTA##*/} >> $LogF/${TopRunName}.txt


# Splitting bamfile in the two 2 amplicon pools, if ref = ampliconref:
if [ "$Reference" == "$AmpliconRef" ]; 
then
	for i in 1 2; do

		# Choose file
		tmpBedFile="${BedFilePool%x.bed}"_pool${i}.bed
		BamInt="${BamFile%.bam}"_intersect${i}.bam
		BamOut="${BamInt%.bam}".readGroupFix.bam
		VarFile=${BamOut%.bam}.vcf

		# Call variants
		gatk --java-options "-Xmx4g" HaplotypeCaller \
		--sample-ploidy 1 \
			-R $Ref_FASTA \
			-I $BamOut \
			-O $VarFile 

		# Filtrerer for bedfil positioner
		gatk --java-options "-Xmx4g" SelectVariants  \
			-R $Ref_FASTA \
			-V $VarFile \
			-L $tmpBedFile \
			-O "${VarFile%.vcf}"_IntFilt.vcf 

	done 

	# Merge vcf files from each pool. Does not merge bam files, as duplicates can occur, because each split can contain some of the same reads:
	java -jar ~/picard.jar MergeVcfs \
	-I "${BamFile%.bam}"_intersect1.readGroupFix_IntFilt.vcf \
	-I "${BamFile%.bam}"_intersect2.readGroupFix_IntFilt.vcf \
	-O "${BamFile%.bam}"_merged.vcf

	bcftools view "${BamFile%.bam}"_merged.vcf | awk '/^#/{print}; !/^#/{if (!uniq[$2]++) print}' > "${BamFile%.bam}"_merged_fix.vcf

	VarFile="${BamFile%.bam}"_merged_fix.vcf

else # Do normal variant call

	# Add tags necessary for GATK because they are not uBAM files 
	# RGLB = Read group library, RGPL = platform, RGPU = platform unit, RGSM = sample name
	java -jar ~/picard.jar AddOrReplaceReadGroups \
	I=$BamFile \
	O="${BamFile%.bam}".readGroupFix.bam \
	RGLB=lib1 \
	RGPL=IonTorrent \
	RGPU=unit1 \
	RGSM=1

	BamOut="${BamFile%.bam}".readGroupFix.bam

	samtools index $BamOut

	VarFile=${BamOut%.bam}.vcf 

	gatk --java-options "-Xmx4g" HaplotypeCaller  \
	--sample-ploidy 1 \
		-R $Ref_FASTA \
		-I $BamOut \
		-O $VarFile 

fi

# Insert mismatch count in summary for all fastq
echo $Reference $(grep -v '^#' $VarFile | wc -l) >> $MMcountFile

# VARIANT FILTERING
# Getting data for R plotting
mkdir -p $currentF/VCFstats
vcfstatF=$currentF/VCFstats
# Getting depth stats
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > $vcfstatF/depth.txt
# Getting QD (quality by depth) stats
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/QD.txt
# Getting Fisher Strand
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;FS=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/FS.txt
# Getting Strand Odds ratio
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;SOR=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/SOR.txt
# Getting root mean square Mapping quality
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;MQ=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/MQ.txt
# Getting MappingQualityRankSumTest (MQRankSum)
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;MQRankSum=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/MQRankSum.txt
# ReadPosRankSumTest (ReadPosRankSum)
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;ReadPosRankSum=\([0-9]*.[0-9]*\);.*$/\1/'  > $vcfstatF/ReadPosRankSum.txt

# Hard filtering variants
gatk --java-options "-Xmx4g" VariantFiltration \
-R $Ref_FASTA \
-V $VarFile \
--filter-expression "QD < 2.0" --filter-name QD \
--filter-expression "DP < 5" --filter-name DP \
--filter-expression "FS > 60.0" --filter-name FS \
--filter-expression "SRQ > 3.0" --filter-name SOR \
--filter-expression "MQ < 50.0" --filter-name MQ \
--filter-expression "MQRankSum < -12.5" --filter-name MQRS \
--filter-expression "ReadPosRankSum < -8.0" --filter-name LoRPRS \
-O "${VarFile%.vcf}"_filtered.vcf

VarFile="${VarFile%.vcf}"_filtered.vcf

# Getting filtered stats
# Getting depth stats
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > $vcfstatF/depth_filt.txt
# Getting QD (quality by depth) stats
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;QD=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/QD_filt.txt
# Getting Fisher Strand
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;FS=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/FS_filt.txt
# Getting Strand Odds ratio
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;SOR=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/SOR_filt.txt
# Getting root mean square Mapping quality
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;MQ=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/MQ_filt.txt
# Getting MappingQualityRankSumTest (MQRankSum)
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;MQRankSum=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/MQRankSum_filt.txt
# ReadPosRankSumTest (ReadPosRankSum)
grep -v "^#" $VarFile | \
cut -f 8 | \
sed 's/^.*;ReadPosRankSum=\([0-9]*.[0-9]*\);.*$/\1/' > $vcfstatF/ReadPosRankSum_filt.txt

# Excluding filtered variants. This is done whether filtered for amplicon positions or not.
gatk --java-options "-Xmx4g" SelectVariants  \
	-R $Ref_FASTA \
	-V $VarFile \
	--exclude-filtered \
	-O "${VarFile%.vcf}"_FiltEx.vcf


# Inserting mismatch count in summary for all alignments, for filtrered variants
echo $Reference $(grep -v '^#' ${VarFile%.vcf}_FiltEx.vcf | wc -l) >> $MMcountFile_filt

sed '/^##source=HaplotypeCaller/d' $VarFile > ${VarFile%.vcf}_headerfix.vcf # removing duplicate "source"
sed '/^##source=HaplotypeCaller/d' ${VarFile%.vcf}_FiltEx.vcf  > ${VarFile%.vcf}_FiltEx_headerfix.vcf # removing duplicate "source"

mv ${VarFile} $currentF/ResultFiles/${FQName}_${Reference}_filt.vcf;
mv ${VarFile}.idx $currentF/ResultFiles/${FQName}_${Reference}.vcf.idx;
rm ${VarFile%.vcf}_FiltEx.vcf ${VarFile%.vcf}_FiltEx.vcf.idx; # Removing files with filtered pos excluded as it was only to count for mismatch file. 

sort -k2 -n $MMcountFile > "${MMcountFile}".sort # Sorting for least mismatches
mv "${MMcountFile}".sort $MMcountFile # Renaming mismatch file

sort -k2 -n $MMcountFile_filt > "${MMcountFile}"_filt.sort # Sorting for least mismatches
mv "${MMcountFile}"_filt.sort $MMcountFile_filt # Renaming mismatch file
