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
        -r|--geno) Reference="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-t|--typetier) typeTier="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -b|--bedfile) BedFileNameX="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-a|--amplicon) AmpliconRef="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

FQName=$(basename ${FastQFile} .fastq)

###################### EDIT ##########################
#Define Folders and params
BedFileName=$BedFileNameX # Regions in bedfile will be soft clipped from bam
BedFilePool=$MainF/References/BedFiles/"${BedFileNameX}"x.bed
#######################################################

# Define locations
QCF=$MainF/QC; 
RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; 
ErrorF=$QCF/Errors; ResultsF=$MainF/Results/$TopRunName	
###################### Find genotypecall ref ##########################

echo -e "Bruger reference" "$Reference" til $FQName

mkdir -p $ResultsF
mkdir -p $ResultsF/$FQName
workD=$ResultsF/$FQName #Overmappe (working directory) med alle bams som hver er mapped til 1 reference

################## BAM CLEANUP & VARIANT DISCOVERY #################
touch $ResultsF/$FQName/MismatchCounts_"${FQName}".txt
touch $ResultsF/$FQName/MismatchCounts_filt_"${FQName}".txt
MMcountFile=$ResultsF/$FQName/MismatchCounts_"${FQName}".txt
MMcountFile_filt=$ResultsF/$FQName/MismatchCounts_filt_"${FQName}".txt # Filtered MM count file

# Laver mappe til hver alignment og aligner og markerer duplicater og indexerer
mkdir -p $workD/"$Reference"
currentF=$workD/"$Reference"
mkdir -p "$currentF"/ResultFiles
Ref_FASTA=$RefF/IndexedRef/"${Reference}"/"${Reference}".fasta # Find reference for picard
echo Aligner nu med reference ${Ref_FASTA##*/}
BamFile=$currentF/${Reference}.bam
# Aligner
bwa mem -t 12 -v 2 -M "$Ref_FASTA" $FastQFile > "${BamFile%bam}"sam # -v = verbosity, 2 for errors og warnings kun
# Sort bamfile
samtools sort "${BamFile%bam}"sam -o "${BamFile%bam}"sort.bam 

java -Xmx28G -jar ~/picard.jar MarkDuplicates -I "${BamFile%bam}"sort.bam -O "${BamFile%bam}"sort.dup.bam \
-M "${BamFile%.bam}"_DupMetrics.txt -REMOVE_DUPLICATES false 2> $ErrorF/markdup_errors_"${Reference}".txt

java -Xmx28G -jar ~/picard.jar MarkDuplicates -I "${BamFile%bam}"sort.bam -O "${BamFile%bam}"sort.dup_rm.bam \
-M "${BamFile%.bam}"_DupMetrics.txt -REMOVE_DUPLICATES true 2> $ErrorF/markdup_errors_"${Reference}".txt

# Rename & index nodupmarked file
samtools index "${BamFile%bam}"sort.bam 
samtools index "${BamFile%bam}"sort.dup.bam
samtools index "${BamFile%bam}"sort.dup_rm.bam

# Her kan angives om dup eller ikke dup skal bruges til resten af kørsel (.dup.)
BamFile="${BamFile%bam}"sort.bam

# Splitter nu bamfil i de 2 amplicon pools, hvis reference er ampliconref:
if [ "$Reference" == "$AmpliconRef" ]; then # || [ $Reference == "${AmpliconRef}"_revised ]
	for i in 1 2; do

		echo Splitting in 2 pools for amplicons

		# Choose file 
		tmpBedFile="${BedFilePool%x.bed}"_pool${i}.bed

		bedtools intersect -a $BamFile -b $tmpBedFile > "${BamFile%.bam}"_intersect${i}.bam
		BamInt="${BamFile%.bam}"_intersect${i}.bam

		clipped_out=$BamInt

		# HaplotypeCaller nødvendigheder
		# Tilføjer tags, nødvendig for GATK, da der ikke arbejdes med uBAM filer 
		# RGLB = Read group library, RGPL = platform, RGPU = platform unit, RGSM = sample name
		java -jar ~/picard.jar AddOrReplaceReadGroups \
		I=$clipped_out \
		O="${clipped_out%.bam}".readGroupFix.bam \
		RGLB=lib1 \
		RGPL=IonTorrent \
		RGPU=unit1 \
		RGSM=1

		BamOut="${clipped_out%.bam}".readGroupFix.bam

		samtools index $BamOut

		VarFile=${BamOut%.bam}.vcf 

		gatk --java-options "-Xmx4g" HaplotypeCaller  \
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

	# MERGER ikke Bam filer fra hver pool, da der så vil opstå duplikater, fordi begge splits kan indeholde nogle af de samme reads 

	# MERGER vcf filer fra hver pool:
	java -jar ~/picard.jar MergeVcfs \
	-I "${BamFile%.bam}"_intersect1.readGroupFix_IntFilt.vcf \
	-I "${BamFile%.bam}"_intersect2.readGroupFix_IntFilt.vcf \
	-O "${BamFile%.bam}"_merged.vcf

	bcftools view "${BamFile%.bam}"_merged.vcf | awk '/^#/{print}; !/^#/{if (!uniq[$2]++) print}' > "${BamFile%.bam}"_merged_fix.vcf

	VarFile="${BamFile%.bam}"_merged_fix.vcf

else # Hvis ref ikke er ampliconref, køres der normal variant calling

	# Tilføjer tags, nødvendig for GATK, da der ikke arbejdes med uBAM filer 
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

# Indsætter mismatch count i samlet fil over alle alignments
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

# Excluderer filtrerede varianter. Dette bliver gjort uanset om der blev filtreret for amplicon position eller ej
gatk --java-options "-Xmx4g" SelectVariants  \
	-R $Ref_FASTA \
	-V $VarFile \
	--exclude-filtered \
	-O "${VarFile%.vcf}"_FiltEx.vcf

#VarFile="${VarFile%.vcf}".filtEx.vcf


# Indsætter mismatch count i samlet fil over alle alignments, for kun de filtrerede varianter
echo $Reference $(grep -v '^#' ${VarFile%.vcf}_FiltEx.vcf | wc -l) >> $MMcountFile_filt

sed '/^##source=HaplotypeCaller/d' $VarFile > ${VarFile%.vcf}_headerfix.vcf # Fjerner duplicate "source"
sed '/^##source=HaplotypeCaller/d' ${VarFile%.vcf}_FiltEx.vcf  > ${VarFile%.vcf}_FiltEx_headerfix.vcf # Fjerner duplicate "source"


# Omdøber slutfil til mere læsbart navn og ligger i slutmappe
# echo cp $currentF/${Reference}.sort.dup.bam $currentF/ResultFiles/${FQName}_${Reference}.sort.bam; 
# echo cp $currentF/${Reference}.sort.dup.bam.bai $currentF/ResultFiles/${FQName}_${Reference}.sort.bam.bai; 
# echo cp ${VarFile%.vcf}_headerfix.vcf $currentF/ResultFiles/${FQName}_${Reference}.vcf;
# echo cp ${VarFile%.vcf}.vcf.idx $currentF/ResultFiles/${FQName}_${Reference}.vcf.idx;
mv $BamOut $currentF/ResultFiles/${FQName}_${Reference}.sort.bam; 
mv ${BamOut}.bai $currentF/ResultFiles/${FQName}_${Reference}.sort.bam.bai; 
mv ${VarFile%.vcf}_headerfix.vcf $currentF/ResultFiles/${FQName}_${Reference}.vcf;
mv ${VarFile%.vcf}.vcf.idx $currentF/ResultFiles/${FQName}_${Reference}.vcf.idx;
rm ${VarFile%.vcf}_FiltEx.vcf ${VarFile%.vcf}_FiltEx.vcf.idx; 

sort -k2 -n $MMcountFile > "${MMcountFile}".sort # Sorting for least mismatches
mv "${MMcountFile}".sort $MMcountFile # Renaming mismatch file

sort -k2 -n $MMcountFile_filt > "${MMcountFile}"_filt.sort # Sorting for least mismatches
mv "${MMcountFile}"_filt.sort $MMcountFile_filt # Renaming mismatch file
