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

#RunName=Pt_81.IonXpress_068_run
#FastQFile=Pt_81.IonXpress_068 #input uden extenstion

###################### EDIT ##########################
#Define Folders and params
MainF=/home/pato/Skrivebord/HPV16_projekt
FastQF=$MainF/FASTQ
MainF=/home/pato/Skrivebord/HPV16_projekt
SNPPOSBedFileName=FASTQfiles_${SuperRunName}_Nuc_change_coords_K02718.1.bed # Regions in bedfile will be soft clipped from bam
SuperRunFolder=$SuperRunName
#######################################################

# Gå til MainF og opret mappe med samme navn som RunName, eks: 
###################### Generate folder structure ####################################
SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder
#####################################################################################
SNPPOSBedFile=$MainF/Annotation_results/$SNPPOSBedFileName
RefList=$(< $RefdF/RefSubtyper_allFastq.txt) # Finder navn på alle referencer

workD=$ResultsF/$RunName #Overmappe (working directory) med alle bams som hver er mapped til 1 reference

for refType in $RefList; do
	
	currentF=$workD/$refType
	Ref_FASTA=$RefF/${refType}/${refType}.fasta # Find reference for picard
	BamFile=$currentF/${refType}.sort.bam

	samtools depth -a $BamFile -b $SNPPOSBedFile > $currentF/SNP_cov.txt

done