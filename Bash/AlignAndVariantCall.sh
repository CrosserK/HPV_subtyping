#!/bin/bash
set -e
set -u
set -o pipefail

# Alle referencer der ligger i $MainF/RefF overmappe (eksl. undermapper) bliver brugt
# Åben script på denne måde:
# HPV_subtyping.sh [Runname] [Fastqname_without_extension]
# Lav indstillinger i EDIT sektion og kør

RunName=$1
FastQFile=$2 #input uden extenstion
SuperRunFolder=$3
MainF=$4
VirRunOut_run=$5
BedFileNameX=$6
AmpliconRef=$7

# HVIS EN SPECIFIK REFERENCE, CTRL+F: "HER DEFINERES REFERENCER"



#RunName=pt_11.IonXpress_040_run
#FastQFile=pt_11.IonXpress_040 #input uden extenstion
#SuperRunFolder=KaroTest_rev2
#MainF=/home/pato/Skrivebord/HPV16_projekt
#Called_subtype=HQ644236.1
#VirRunOut_run=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/KaroTest/pt_40.IonXpress_065_run
#BedFileNameX=Revised_IAD209923_226_Designed_compl
#AmpliconRef=K02718.1



#######TEST#####
#echo $RunName $FastQFile $SuperRunFolder $MainF $VirRunOut_run $BedFileNameX
#echo og Called_subtype er $Called_subtype
################






###################### EDIT ##########################
#Define Folders and params
BedFileName=$( echo  ${BedFileNameX}) # Regions in bedfile will be soft clipped from bam
BedFilePool=$MainF/References/BedFiles/${BedFileNameX}x.bed
#######################################################

MainF=/home/pato/Skrivebord/HPV16_projekt
SuperRunName=$SuperRunFolder
VirSupOut=$MainF/VirStrain_run/$SuperRunName
VirRunOut_run=$VirSupOut/$RunName

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder

# Fetching Bed file
BedFile=$RefF/$BedFileName
###################### Create directory for each reference ##########################
# Finder navne for hver ref i References overmappe (eksl. undermapper)
find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_${RunName}.txt
RefListOrigin=$(< $RefdF/RefSubtyper_${RunName}.txt) # Finder navn på alle referencer i mappe

# HER DEFINERES REFERENCER
# RefListCalls=K02718.1_revised
RefListCalls=$(< $VirRunOut_run/SubtypeCall.txt)

echo Har kald: $RefListCalls slut

# Clearer RefList for referencer som ikke er i RefList mappe (eks fordi der ikke kunne kaldes op til 3 mulige typer og der derfor er blevet trukket forkert string ud)
touch $VirRunOut_run/Revised_SubTypeCalls.txt
RevRefs=$VirRunOut_run/Revised_SubTypeCalls.txt


############## Clear unviable subtype calls #################
for refType in $RefListOrigin; do

for refCall in $RefListCalls; do

if [ $refType == $refCall ]; then

Check=$(grep "$refCall" $VirRunOut_run/SubtypeCall.txt)

if [ ${#Check} -gt 0 ]; then
	
	echo $refCall >> $RevRefs
fi

fi
done
done
# unset Called_subtype
##############################################################

# Tager nu revideret subtype fil:
RefList=$(< $RevRefs)
export RevRefCalls=$(< $RevRefs) # Giver revideret refs til main script
#### TEST ####
#echo Til $FastQFile bruger fil $RefdF/RefSubtyper_${RunName}.txt og har revideret til "$RefList"

######################################################################################

############### FASTQ preprocessing and filtering (DATA CLEANUP) ##################### 
# Se read længder og antal reads med awk
#cat $SeqF/${FastQFile}.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c  # Read lenghts, 1. kolonne er antal reads, 2. kolonne er længde
# echo $(cat $SeqF/${FastQFile}.fastq | wc -l)/4 | bc # Number of reads
######################################################################################

####################### FASTQ QC & Filter #############################
# # -O Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
#fastqc -o $QualF $SeqF/${FastQFile}.fastq 
# cutadapt -q $QualTrim -m $cutadaptMinsize -M $cutadaptMaxsize $SeqF/${FastQFile}.fastq -o $SeqF/${FastQFile}_filt.fastq
#fastqc -o $QualF $SeqF/${FastQFile}_filt.fastq 
# firefox $QualF $SeqF/${FastQFile}_filt_fastqc.html &
##################### DEFINING RUNFOLDER ######################

mkdir -p $ResultsF
mkdir -p $ResultsF/$RunName
workD=$ResultsF/$RunName #Overmappe (working directory) med alle bams som hver er mapped til 1 reference

################## BAM CLEANUP & VARIANT DISCOVERY #################
touch $ResultsF/$RunName/MismatchCounts_${RunName}.txt
touch $ResultsF/$RunName/MismatchCounts_filt_${RunName}.txt
MMcountFile=$ResultsF/$RunName/MismatchCounts_${RunName}.txt
MMcountFile_filt=$ResultsF/$RunName/MismatchCounts_filt_${RunName}.txt # Filtered MM count file

for refType in $RefList; do 

	####TEST
	#refType=AF536179.1
	# Laver mappe til hver alignment og aligner og markerer duplicater og indexerer
	mkdir -p $workD/$refType
	currentF=$workD/$refType
	mkdir -p $currentF/ResultFiles
	Ref_FASTA=$RefF/${refType}/${refType}.fasta # Find reference for picard
	echo reference er nu $Ref_FASTA
	BamFile=$currentF/${refType}.bam
	# Aligner
	bwa mem -t 12 -v 2 $Ref_FASTA $SeqF/${FastQFile}_filt.fastq > ${BamFile%bam}sam # -v = verbosity, 2 for errors og warnings kun

	# Sort bamfile
	samtools sort ${BamFile%bam}sam -o ${BamFile%bam}sort.bam 

	java -Xmx28G -jar ~/picard.jar MarkDuplicates -I ${BamFile%bam}sort.bam -O ${BamFile%bam}sort.dup.bam \
	-M $currentF/${refType}_DupMetrics.txt -REMOVE_DUPLICATES false 2> $ErrorF/markdup_errors_${refType}.txt

	java -Xmx28G -jar ~/picard.jar MarkDuplicates -I ${BamFile%bam}sort.bam -O ${BamFile%bam}sort.dup_rm.bam \
	-M $currentF/${refType}_DupMetrics.txt -REMOVE_DUPLICATES true 2> $ErrorF/markdup_errors_${refType}.txt

	# Rename & index nodupmarked file
	samtools index ${BamFile%bam}sort.bam 
	samtools index ${BamFile%bam}sort.dup.bam
	samtools index ${BamFile%bam}sort.dup_rm.bam

	BamFile=${BamFile%bam}sort.dup.bam

	# Splitter nu bamfil i de 2 amplicon pools, hvis reference er ampliconref:
	if [ $refType == $AmpliconRef ] || [ $refType == ${AmpliconRef}_revised ]; then
		for i in 1 2; do

			# Choose file 
			tmpBedFile=${BedFilePool%x.bed}${i}.bed

			bedtools intersect -a $BamFile -b $tmpBedFile > ${BamFile%.bam}_intersect${i}.bam
			BamInt=${BamFile%.bam}_intersect${i}.bam

			clipped_out=$BamInt

			# HaplotypeCaller nødvendigheder
			# Tilføjer tags, nødvendig for GATK, da der ikke arbejdes med uBAM filer 
			# RGLB = Read group library, RGPL = platform, RGPU = platform unit, RGSM = sample name
			java -jar picard.jar AddOrReplaceReadGroups \
			I=$clipped_out \
			O=${clipped_out%.bam}.readGroupFix.bam \
			RGLB=lib1 \
			RGPL=IonTorrent \
			RGPU=unit1 \
			RGSM=1

			BamOut=${clipped_out%.bam}.readGroupFix.bam

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
			   -O ${VarFile%.vcf}_IntFilt.vcf 

		done 

		# MERGER ikke Bam filer fra hver pool, da der så vil opstå duplikater, fordi begge splits kan indeholde nogle af de samme reads
	
		# MERGER vcf filer fra hver pool:
	   	java -jar picard.jar MergeVcfs \
	   	-I ${BamFile%.bam}_intersect1.readGroupFix_IntFilt.vcf \
	   	-I ${BamFile%.bam}_intersect2.readGroupFix_IntFilt.vcf \
	   	-O ${BamFile%.bam}_merged.vcf

		bcftools view ${BamFile%.bam}_merged.vcf | awk '/^#/{print}; !/^#/{if (!uniq[$2]++) print}' > ${BamFile%.bam}_merged_fix.vcf

		VarFile=${BamFile%.bam}_merged_fix.vcf
	
	else # Hvis ref ikke er ampliconref, køres der normal variant calling

		#### TEST
		echo splitter ikke i 2 pools

		# HaplotypeCaller nødvendigheder
		# Tilføjer tags, nødvendig for GATK, da der ikke arbejdes med uBAM filer 
		# RGLB = Read group library, RGPL = platform, RGPU = platform unit, RGSM = sample name
		java -jar picard.jar AddOrReplaceReadGroups \
		I=$BamFile \
		O=${BamFile%.bam}.readGroupFix.bam \
		RGLB=lib1 \
		RGPL=IonTorrent \
		RGPU=unit1 \
		RGSM=1

		BamOut=${BamFile%.bam}.readGroupFix.bam

		samtools index $BamOut

		VarFile=${BamOut%.bam}.vcf 

		gatk --java-options "-Xmx4g" HaplotypeCaller  \
		--sample-ploidy 1 \
		   -R $Ref_FASTA \
		   -I $BamOut \
		   -O $VarFile 

		## Filtrerer for bedfil positioner
		#gatk --java-options "-Xmx4g" SelectVariants  \
		#   -R $Ref_FASTA \
		#   -V $VarFile \
		#   -L $tmpBedFile \
		#   -O ${VarFile%.vcf}_IntFilt.vcf 

		#VarFile=${VarFile%.vcf}_IntFilt.vcf 


	fi
	echo nu skal der laves $currentF/VCFstats

	
	# Indsætter mismatch count i samlet fil over alle alignments
	echo $refType $(grep -v '^#' $VarFile | wc -l) >> $MMcountFile

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
   --filter-expression "MQ < 40.0" --filter-name MQ \
   --filter-expression "MQRankSum < -12.5" --filter-name MQRS \
   --filter-expression "ReadPosRankSum < -8.0" --filter-name LoRPRS \
   -O ${VarFile%.vcf}_filtered.vcf

   VarFile=${VarFile%.vcf}_filtered.vcf

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


    #VarFile=${VarFile%.vcf}.filtEx.vcf


	# Indsætter mismatch count i samlet fil over alle alignments
	echo $refType $(grep -v '^#' $VarFile | wc -l) >> $MMcountFile_filt

	sed '/^##source=HaplotypeCaller/d' $VarFile > ${VarFile%.vcf}_headerfix.vcf

	# Omdøber slutfil til mere læsbart navn og ligger i slutmappe
	mv $currentF/${refType}.sort.dup.bam $currentF/ResultFiles/${FastQFile}_${refType}.sort.bam 
	mv $currentF/${refType}.sort.dup.bam.bai $currentF/ResultFiles/${FastQFile}_${refType}.sort.bam.bai 
	mv ${VarFile%.vcf}_headerfix.vcf $currentF/ResultFiles/${FastQFile}_${refType}.vcf
	mv ${VarFile%.vcf}.vcf.idx $currentF/ResultFiles/${FastQFile}_${refType}.vcf.idx



done

unset RefList

sort -k2 -n $MMcountFile > ${MMcountFile}.sort # Sorting for least mismatches
mv ${MMcountFile}.sort $MMcountFile # Renaming mismatch file

sort -k2 -n $MMcountFile_filt > ${MMcountFile}_filt.sort # Sorting for least mismatches
mv ${MMcountFile}_filt.sort $MMcountFile_filt # Renaming mismatch file

# Rydder op
# Sikrer at RefF er korrekt angivet, så der ikke slettes for meget
# (i tilfælde af at den er blevet slettet længere oppe, eller RunName variabel er tom)
# -gt står for greater than
RefF=$MainF/References/$RunName 
#if [ ${#RunName} -gt 0 ]; then
#	rm -r $RefF
#fi
rm $RefdF/RefSubtyper_${RunName}.txt
# Fjerner filtreret fastq, da den ikke skal bruges længere
#rm $SeqF/${FastQFile}_filt.fastq




