#!/bin/bash
set -e
set -u
set -o pipefail


# Generate indices for HPV_subtyping.sh
# Alle referencer der ligger i $MainF/RefF overmappe (eksl. undermapper) bliver brugt
# Åben script på denne måde:
# Idx_for_HPV_subtyping.sh [SuperRunName] 
# Lav indstillinger i EDIT sektion og kør

SuperRunName=$1

###################### EDIT ##########################
#Define Folders and params
MainF=/home/pato/Skrivebord/HPV16_projekt
#######################################################

# Gå til MainF og opret mappe med samme navn som SuperRunName, eks: 
###################### Generate folder structure ####################################
mkdir -p $MainF/{Aligned/{Samfiles,BamFiles/unsorted},References,Results/$SuperRunName,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}}
# -p option enables returning no error if folders exist. They will not be overwritten either way. 
SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results
#####################################################################################
# Fetching Bed file
BedFile=$RefF/$BedFileName
###################### Create directory for each reference ##########################
# Finder navne for hver ref i References overmappe (eksl. undermapper)
find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_${SuperRunName}.txt
RefList=$(< $RefdF/RefSubtyper_${SuperRunName}.txt) # Finder navn på alle referencer



 # Putter hver reference i en undermappe for sig, så de er klar til kørsel
 # genererer også dict fil og samtools faidx fil til HaplotypeCaller
 # Genererer index til bwa mem
for refType in $RefList; do
	
	mkdir -p $RefF/$SuperRunName/$refType
	cp $RefF/${refType}.fasta $RefF/$SuperRunName/${refType}/${refType}.fasta

	Ref_FASTA=$RefF/$SuperRunName/${refType}/${refType}.fasta

	# Create sequence dictionary for gatk haplotypecaller
	java -jar picard.jar CreateSequenceDictionary \
	      R=$Ref_FASTA \
	      O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt

	# Create index for bwa mem
	bwa index $Ref_FASTA

	# Index with samtools faidx
	samtools faidx $Ref_FASTA

done
