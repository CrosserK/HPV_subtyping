#!/bin/bash
set -e
set -u
set -o pipefail

# Åben script på denne måde:
# Specific_site_cov.sh [Runname] [Fastqname_without_extension]
# Lav indstillinger i EDIT sektion og kør

RunName=$1
FastQFile=$2 #input uden extenstion
SuperRunName=$3

#RunName=pt_40.IonXpress_065_run
#FastQFile=pt_40.IonXpress_065 #input uden extenstion
#SuperRunName=Karoline_run_test2


RefListFile=/home/pato/Skrivebord/HPV16_projekt/ReferenceDetails/RefSubtyper_latest.txt


#Define Folders and params
MainF=/home/pato/Skrivebord/HPV16_projekt
FastQF=$MainF/FASTQ


SuperRunFolder=$SuperRunName
# Gå til MainF og opret mappe med samme navn som RunName, eks: 
###################### Generate folder structure ####################################
SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder
#####################################################################################
#RefList=$(< $RefdF/${RefListFile}.txt) # Finder navn på alle referencer
RefList=$(< $RefListFile) 

workD=$ResultsF/$RunName #Overmappe (working directory) med alle bams som hver er mapped til 1 reference

for refType in $RefList; do
	if [ -d $ResultsF/$RunName/$refType ]; then
	###################### EDIT ##########################
	SNPPOSBedFileName=$(echo FASTQfiles_${SuperRunName}_Nuc_change_coords_${refType}.bed) # Regioner der skal tjekkes
	#######################################################
	SNPPOSBedFile=$MainF/Annotation_results/$SNPPOSBedFileName
	
	currentF=$workD/$refType/ResultFiles
	Ref_FASTA=$RefF/IndexedRef/${refType}/${refType}.fasta # Find reference for picard
	BamFile=$currentF/${FastQFile}_${refType}.sort.bam

	samtools depth -a $BamFile -b $SNPPOSBedFile > $currentF/SNP_cov.txt
	fi
done