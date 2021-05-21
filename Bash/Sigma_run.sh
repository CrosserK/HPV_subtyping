#!/bin/bash
set -e
set -u
set -o pipefail

# Åben script på denne måde:
# Sigma_run.sh [Runname] [Fastqname_without_extension]
# Lav indstillinger i EDIT sektion og kør

###################### EDIT ##########################
RunName=$1
FastQFile=$2 #input uden .fastq

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
#######################################################


###################### Generate folder structure ####################################
mkdir -p $MainF/{Aligned/{Samfiles,BamFiles/unsorted},References,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics}} 
# -p option enables returning no error if folders exist. They will not be overwritten either way. 

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth
FlagF=$AnaF/Flagstats; DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References
#####################################################################################


###################### Create directory for each reference ##########################
# Finder navne for hver ref i References overmappe (eksl. undermapper)
find $RefF/ -maxdepth 1  -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_${RunName}.txt


RefList=$(< $RefdF/RefSubtyper_${RunName}.txt)  
 # Putter hver reference i en undermappe for sig, så de er klar til Sigma, som også kan indexere dem
for refType in $RefList; do
	mkdir -p $RefF/$RunName/$refType
	cp $RefF/${refType}.fasta $RefF/$RunName/${refType}/${refType}.fasta
done

# Bruger nu genereret undermappe med referencer
RefF=$MainF/References/$RunName

######################################################################################



#################### FASTQ preprocessing and filtering (DATA CLEANUP) ################

# Se read længder og antal reads med awk
#cat $SeqF/${FastQFile}.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c  # Read lenghts, 1. kolonne er antal reads, 2. kolonne er længde
# echo $(cat $SeqF/${FastQFile}.fastq | wc -l)/4 | bc # Number of reads

#fastqc -o $QualF $SeqF/${FastQFile}.fastq 
cutadapt -O 3 -q 20 -m 75 -M 150 $SeqF/${FastQFile}.fastq -o $SeqF/${FastQFile}_filt.fastq

######################################################################################



####################### FASTQ QC #############################

# -m er minumum længde, -M er maksimum længde
# # Fjerner low qual baser. -q betyder at den fjerner baser under den kvalitet fra 3' enden. 
# Hvis der angives et tal mere separeret af komma (-q 20,20) vil den også trimme fra 5' enden
# # -m [N] minimum lenght. -M [N] maximum lenght
# # -O Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3
# fastqc -o $QualF $SeqF/${FastQFile}_trim.fastq 
# seqkit seq $SeqF/${FastQFile}_trim.fastq -m 75 -M 150 > $SeqF/${FastQFile}_filt.fastq # Size filtering
# fastqc -o $QualF $SeqF/${FastQFile}_filt.fastq 
# firefox $QualF $SeqF/${FastQFile}_filt_fastqc.html &

#mv $SeqF/${FastQFile}.fastq $SeqF/${FastQFile}_untouched.fastq
mv $SeqF/${FastQFile}_filt.fastq $SeqF/${FastQFile}.fastq

#fastqc -o $QualF $SeqF/${FastQFile}.fastq
#firefox $QualF $SeqF/${FastQFile}_fastqc.html &

###############################################################

##################### DEFINING RUNFOLDER ######################

mkdir -p $SigmaF
mkdir -p $SigmaF/$RunName
BamF=$SigmaF/$RunName #Overmappe med alle bams som hver er mapped til 1 reference

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
sigma-index-genomes -p 12 -c $SigmaF/sigma_config_${RunName}.cfg -w $BamF
# Bowtie antal threads er defineret i configfil, hvis yderligere tilføjes her, multipleres antal threads forsøgt benyttet
sigma-align-reads -c $SigmaF/sigma_config_${RunName}.cfg -w $BamF

#######################################################

################## BAM CLEANUP & VARIANT DISCOVERY #################

SigmaBamF=$SigmaF/$RunName/sigma_alignments_output

touch $SigmaF/$RunName/MismatchCounts.txt
MMcountFile=$SigmaF/$RunName/MismatchCounts.txt
# Remove duplicate reads and rename to original name, store original file with undup in name. Sort and index bam file

for refType in $RefList; do 
	currentF=$SigmaBamF/$refType
	

	
	samtools sort $currentF/${refType}.align.bam -o $currentF/${refType}.align.sort.bam

	java -Xmx28G -jar ~/picard.jar MarkDuplicates -I $currentF/${refType}.align.sort.bam -O $currentF/${refType}.align_dup.bam -M $currentF/${refType}_DupMetrics.txt -REMOVE_DUPLICATES false
    mv $currentF/${refType}.align.bam $currentF/${refType}.align.undup.bam
    mv $currentF/${refType}.align_dup.bam $currentF/${refType}.align.bam

    # Indel realignment

    # Base recalibration

    samtools index $currentF/${refType}.align.bam 
    Ref_FASTA=$RefF/${refType}/${refType}.fasta
    VarFile=$currentF/${refType}.align.vcf
	freebayes -p 1 -f $Ref_FASTA $currentF/${refType}.align.bam > $VarFile # -C for at bestemme hvor mange reads skal support en variant. -p for ploidy. Freebayes inkluderer leftalignment
	# Indsætter mismatch count i fil
	echo $refType $(grep -v '^#' $VarFile | wc -l) >> $MMcountFile
	sort -k2 -n $MMcountFile > ${MMcountFile}.sort
	mv ${MMcountFile}.sort $MMcountFile

	# Variant recalibration


done

# Denne kører både build-model og solve-model moduler. 2 fastq på 10 HPV subtyper = 383.835s
# Nogle af filerne kan godt ende i home/pato, selvom working directory er defineret. Søg efter "sigma_out"
sigma -t 12 -c $SigmaF/sigma_config_${RunName}.cfg -w $BamF

mv sigma_out.html $SigmaF/$RunName/sigma_out.html
mv sigma_out.ipopt.txt $SigmaF/$RunName/sigma_out.ipopt.txt
mv sigma_out.gvector.txt $SigmaF/$RunName/sigma_out.gvector.txt ;








## Ekstra analyse
#sigma-bootstrap -t 12 –p 12 -w $BamF

#sigma-jackknife -t 12 –p 12 -w $BamF

# Opdater variant calling i config fil
#sigma-target-reads -t 12 -p 12 -w $BamF