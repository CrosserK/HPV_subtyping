#!/bin/bash
set -e
set -u
set -o pipefail

# Sigma config_init.cfg skal ligge i SigmaF mappe
# Fastq fil path og ref path skal defineres ved sed linje under edit, selve navnet bliver hentet fra input $2
# Alle referencer der ligger i RefF overmappe (eksl. undermapper) bliver brugt
# Åben script på denne måde:
# Sigma_run.sh [Runname] [Fastqname_without_extension]
# Lav indstillinger i EDIT sektion og kør

RunName=$1
FastQFile=$2 #input uden extenstion

###################### EDIT ##########################
#Define Folders
MainF=/home/pato/Skrivebord/HPV16_projekt
SigmaF=$MainF/Sigma_run

# Gå til Sigma_run og opret mappe med samme navn som RunName, eks: 
mkdir $SigmaF/$RunName
# kopier derefter configfil fra tidligere kørsel i mappen ovenover og giv navn som her:
cp $SigmaF/sigma_config_init.cfg $SigmaF/sigma_config_${RunName}.cfg
# og edit den med fastq fil (VIGTIGT AT DEN HAR SIN EXTENSION .fastq på) og reference mappe (skal inkludere run name eks:
# ~/HPV16_projekt/References/run04
sed -i "29s/.*/Reference_Genome_Directory=\/home\/pato\/Skrivebord\/HPV16_projekt\/References\/$RunName/" $SigmaF/sigma_config_${RunName}.cfg 
# Referencer der ligger frit i References mappe (eksl. undermapper) bliver brugt og lagt i hver deres undermappe for at Sigma kan benytte dem
sed -i "41s/.*/Single_End_Reads=\/home\/pato\/Skrivebord\/HPV16_projekt\/FASTQ\/${FastQFile}.fastq/" $SigmaF/sigma_config_${RunName}.cfg
cutadaptMinsize=50
cutadaptMaxsize=320
#######################################################


###################### Generate folder structure ####################################
mkdir -p $MainF/{Aligned/{Samfiles,BamFiles/unsorted},References,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}} 
# -p option enables returning no error if folders exist. They will not be overwritten either way. 

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References ; ErrorF=$AnaF/Errors

#####################################################################################
# Fetching Bed file
BedFile=$RefF/IAD209923_226_Designed.bed

###################### Create directory for each reference ##########################
# Finder navne for hver ref i References overmappe (eksl. undermapper)
find $RefF/ -maxdepth 1  -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_${RunName}.txt


RefList=$(< $RefdF/RefSubtyper_${RunName}.txt) # Finder navn på alle referencer

 # Putter hver reference i en undermappe for sig, så de er klar til Sigma, som også kan indexere dem til bowtie2, 
 # genererer også dict fil og samtools faidx fil til HaplotypeCaller
 # Genererer index til bwa mem
for refType in $RefList; do
	
	mkdir -p $RefF/$RunName/$refType
	cp $RefF/${refType}.fasta $RefF/$RunName/${refType}/${refType}.fasta

	Ref_FASTA=$RefF/$RunName/${refType}/${refType}.fasta

	# Create sequence dictionary for gatk haplotypecaller
	java -jar picard.jar CreateSequenceDictionary \
	      R=$Ref_FASTA \
	      O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt

	# Create index for bwa mem
	bwa index $Ref_FASTA

	# Index with samtools faidx
	samtools faidx $Ref_FASTA

done

# Bruger nu genereret undermappe med referencer
RefF=$MainF/References/$RunName
######################################################################################



#################### FASTQ preprocessing and filtering (DATA CLEANUP) ################
# Se read længder og antal reads med awk
#cat $SeqF/${FastQFile}.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c  # Read lenghts, 1. kolonne er antal reads, 2. kolonne er længde
# echo $(cat $SeqF/${FastQFile}.fastq | wc -l)/4 | bc # Number of reads
######################################################################################



####################### FASTQ QC & Filter #############################

# # -O Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
fastqc -o $QualF $SeqF/${FastQFile}.fastq 
cutadapt -q 20 -m $cutadaptMinsize -M $cutadaptMaxsize $SeqF/${FastQFile}.fastq -o $SeqF/${FastQFile}_filt.fastq
fastqc -o $QualF $SeqF/${FastQFile}_filt.fastq 
# firefox $QualF $SeqF/${FastQFile}_filt_fastqc.html &

#mv $SeqF/${FastQFile}.fastq $SeqF/${FastQFile}_untouched.fastq
# mv $SeqF/${FastQFile}_filt.fastq $SeqF/${FastQFile}.fastq # Omdøber filteret fil til start navn, original fil bliver slettet
###############################################################

##################### DEFINING RUNFOLDER ######################

mkdir -p $SigmaF
mkdir -p $SigmaF/$RunName
workD=$SigmaF/$RunName #Overmappe (working directory) med alle bams som hver er mapped til 1 reference

# # Get names of all references ind subdirectories
#touch $MainF/ReferenceDetails/RefSubtyper_${RunName}_temp.txt 
#RefSubtyper=RefSubtyper_${RunName} 

#Hvis referencer allerede ligger i undermapper:
# for refType in $(($RefdF/HPV_subtypeList.txt)) ; do
#     echo "$refType" >> $MainF/{$RefSubtyper}_${RunName}_temp.txt
# done

#tr -d \/ < $MainF/{$RefSubtyper}_${RunName}_temp.txt > $MainF/{$RefSubtyper}_${RunName}_temp1.txt
#rm {$RefSubtyper}_${RunName}_temp.txt

#sed -n '/References//,$p' $MainF/{$RefSubtyper}_${RunName}_temp1.txt > $MainF/{$RefSubtyper}_${RunName}.txt
#Skal være alle som ligger i mappernes

##############################################################



################### ALIGNMENT #########################

# Remove duplicate, identifying and masking human reads, trim low qual bases 
# Sigma index genomes for bowtie2
#sigma-index-genomes [options] -c <config file path> -w <working directory>
# sigma-index-genomes -p 12 -c $SigmaF/sigma_config_${RunName}.cfg -w $workD
# Bowtie antal threads er defineret i configfil, hvis yderligere tilføjes her, multipleres antal threads forsøgt benyttet
# sigma-align-reads -c $SigmaF/sigma_config_${RunName}.cfg -w $workD

#######################################################

################## BAM CLEANUP & VARIANT DISCOVERY #################

SigmaBamF=$SigmaF/$RunName/sigma_alignments_output
MMcountFile=$SigmaF/$RunName/MismatchCounts_${RunName}.txt
MMcountFile_filt=$SigmaF/$RunName/MismatchCounts_filt_${RunName}.txt # Filtered MM count file

# For hver bamfil alignet til unik reference, fortag nu:

## python 'VirStrain_build.py' -i '/home/pato/Skrivebord/HPV16_projekt/References_andre/15_HPV_sublineages.mafft' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_15_virstrain'

for refType in $RefList; do 

	# Laver mappe til hver alignment
	mkdir -p $SigmaBamF/$refType
	currentF=$SigmaBamF/$refType
	Ref_FASTA=$RefF/${refType}/${refType}.fasta # Find reference for picard
	BamFile=$currentF/${refType}.align.bam
	# Aligner
	bwa mem -t 12 -v 2 $Ref_FASTA $SeqF/${FastQFile}_filt.fastq > ${BamFile%bam}sam # -v = verbosity, 2 for errors og warnings kun

	# Sort bamfile
	samtools sort ${BamFile%bam}sam -o ${BamFile%bam}sort.bam 

	# Lft align
	#gatk LeftAlignIndels \
   #-R $Ref_FASTA \
   #-I ${BamFile%bam}sort.bam \
   #-O ${BamFile%bam}sort.leftAl.bam

   	#mv ${BamFile%bam}sort.bam ${BamFile%bam}sort.noLeftAl.bam
   	#samtools index ${BamFile%bam}sort.noLeftAl.bam
   	#mv ${BamFile%bam}sort.leftAl.bam ${BamFile%bam}sort.bam

	java -Xmx28G -jar ~/picard.jar MarkDuplicates -I ${BamFile%bam}sort.bam -O ${BamFile%bam}sort.dup.bam \
	-M $currentF/${refType}_DupMetrics.txt -REMOVE_DUPLICATES false 2> $ErrorF/markdup_errors_${refType}.txt

	java -Xmx28G -jar ~/picard.jar MarkDuplicates -I ${BamFile%bam}sort.bam -O ${BamFile%bam}sort.dup_rm.bam \
	-M $currentF/${refType}_DupMetrics.txt -REMOVE_DUPLICATES true 2> $ErrorF/markdup_errors_${refType}.txt

    # Rename & index nodupmarked file
    samtools index ${BamFile%bam}sort.bam 
    samtools index ${BamFile%bam}sort.dup.bam
    samtools index ${BamFile%bam}sort.dup_rm.bam

    # HaplotypeCaller nødvendigheder
	# Tilføjer tags, nødvendig for GATK, da der ikke arbejdes med uBAM filer (sigma skal bruge fastq filer)
	java -jar picard.jar AddOrReplaceReadGroups \
	I=${BamFile%bam}sort.dup.bam \
	O=${BamFile%bam}sort.dup.readGroupFix.bam \
	RGLB=lib1 \
	RGPL=IonTorrent \
	RGPU=unit1 \
	RGSM=1
	# RGLB = Read group library, RGPL = platform, RGPU = platform unit, RGSM = sample name

	# Index bam
    samtools index ${BamFile%bam}sort.dup.readGroupFix.bam

	# Intersecting bam with bed (removing everything not in bed)
	bedtools intersect -a ${BamFile%bam}sort.dup.readGroupFix.bam -b $BedFile > ${BamFile%bam}sort.dup.readGroupFix.isec.bam

	samtools index ${BamFile%bam}sort.dup.readGroupFix.isec.bam

	VarFile=$currentF/${refType}.align.vcf # Define vcf name

	gatk --java-options "-Xmx4g" HaplotypeCaller  \
	--sample-ploidy 1 \
	   -R $Ref_FASTA \
	   -I ${BamFile%bam}sort.dup.readGroupFix.isec.bam\
	   -O $VarFile 

	freebayes -p 1 -f $Ref_FASTA ${BamFile%bam}sort.dup.readGroupFix.isec.bam > ${VarFile%vcf}freeb.vcf # -C for at bestemme hvor mange reads skal support en variant. -p for ploidy. Freebayes inkluderer leftalignment
	mv ${VarFile%vcf}freeb.vcf $VarFile

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

	# Hard filtering variants. Note that FS is not important for non-poliploid organisms
	gatk --java-options "-Xmx4g" VariantFiltration \
   -R $Ref_FASTA \
   -V $VarFile \
   --filter-expression "QD < 2.0" --filter-name QD \
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

	# Removing filtered sites
	gatk --java-options "-Xmx4g" SelectVariants \
     -R $Ref_FASTA \
     -V $VarFile \
     --exclude-filtered \
     -O ${VarFile%.vcf}.filtEx.vcf # Filtered exluded from vcf

    VarFile=${VarFile%.vcf}.filtEx.vcf

	# Indsætter mismatch count i samlet fil over alle alignments
	echo $refType $(grep -v '^#' $VarFile | wc -l) >> $MMcountFile_filt

done

sort -k2 -n $MMcountFile > ${MMcountFile}.sort # Sorting for least mismatches
mv ${MMcountFile}.sort $MMcountFile # Renaming mismatch file

sort -k2 -n $MMcountFile_filt > ${MMcountFile}_filt.sort # Sorting for least mismatches
mv ${MMcountFile}_filt.sort $MMcountFile_filt # Renaming mismatch file

# View bam stats command line interface
# java -Xmx4g -jar ./BAMStats-1.25/BAMStats-1.25.jar -i $BamFile
# View bam stats GUI
# java -Xmx4g -jar ./BAMStats-1.25/BAMStats-GUI-1.25.jar -i $BamFile


# Denne kører både build-model og solve-model moduler. 2 fastq på 10 HPV subtyper = 383.835s
# Nogle af filerne kan godt ende i home/pato, selvom working directory er defineret. Søg efter "sigma_out"
# !!! sigma -t 12 -c $SigmaF/sigma_config_${RunName}.cfg -w $workD

# !!! mv sigma_out.html $SigmaF/$RunName/sigma_out.html
# !!! mv sigma_out.ipopt.txt $SigmaF/$RunName/sigma_out.ipopt.txt
# !!! mv sigma_out.gvector.txt $SigmaF/$RunName/sigma_out.gvector.txt ;

# Rydder op
rm $SigmaF/sigma_config_${RunName}.cfg
RefF=$MainF/References/$RunName # Sikrer at RefF er korrekt angivet, så der ikke slettes for meget
# (i tilfælde af at den er blevet slettet længere oppe)
rm -r $RefF




## Ekstra analyse
#sigma-bootstrap -t 12 –p 12 -w $workD

#sigma-jackknife -t 12 –p 12 -w $workD

# Opdater variant calling i config fil
#sigma-target-reads -t 12 -p 12 -w $workD