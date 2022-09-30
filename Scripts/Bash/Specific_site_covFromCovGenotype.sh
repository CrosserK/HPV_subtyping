#!/bin/bash
set -e
set -u
set -o pipefail

# Lav indstillinger i EDIT sektion og kør
# Åben script på denne måde:
# Specific_site_cov.sh [Runname] [Fastqname_without_extension]

MainF=$1
RunName=$2
FastQFile=$3 #input uden extenstion
SuperRunName=$4

#MainF=/home/pato/Skrivebord/HPV_subtyping
#RunName=pt_138.IonXpress_007_run
#FastQFile=pt_138.IonXpress_007 #input uden extenstion
#SuperRunName=Karoline_run_covGenTest_0859_07072021

#Define Folders and params
FastQF=$MainF/FASTQ


SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunName
#####################################################################################

# Tjekker om der er splittede bam filer
if [ -f "$MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt" ]; then
	RefList=$(< $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt) # Finder navn på alle referencer
else
	RefList=$(< $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt) # Finder navn på alle referencer
fi

workD=$ResultsF/$RunName #Overmappe (working directory) med alle bams som hver er mapped til 1 reference

for refType in $RefList; do

	SNPPOSBedFileName=$(echo FASTQfiles_${SuperRunName}_Nuc_change_coords_${refType}.bed) # Regioner der skal tjekkes
	SNPPOSBedFile=$MainF/Annotation_results/$SNPPOSBedFileName
	currentF=$workD/$refType/ResultFiles
	Ref_FASTA=$RefF/IndexedRef/${refType}/${refType}.fasta # Find reference for picard
	BamFile=$currentF/${FastQFile}_${refType}.sort.bam

	# Tjekker om der findes en SNPPosBedFil, ellers er det fordi der ikke var nogen varianter i den HPVtype
	if [ -f $SNPPOSBedFile ]; then
	samtools depth -a $BamFile -b $SNPPOSBedFile > $currentF/SNP_cov.txt
	fi
done