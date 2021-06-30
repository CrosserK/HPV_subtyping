#!/bin/bash
set -e
set -u
#set -o pipefail

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
customRef=true

# TEST
#RunName=pt_7.IonXpress_038_run
#FastQFile=pt_7.IonXpress_038 #input uden extenstion
#SuperRunFolder=Karoline_run_test
#MainF=/home/pato/Skrivebord/HPV16_projekt
#VirRunOut_run=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/"{SuperRunFolder}"/"${RunName}"
#BedFileNameX=IAD209923_226_Designed_compl
#AmpliconRef=K02718.1

#######TEST#####
#echo $RunName $FastQFile $SuperRunFolder $MainF $VirRunOut_run $BedFileNameX
#echo og Called_subtype er $Called_subtype
################



###################### EDIT ##########################
#Define Folders and params
BedFileName=$BedFileNameX # Regions in bedfile will be soft clipped from bam
BedFilePool=$MainF/References/BedFiles/"${BedFileNameX}"x.bed
#######################################################

MainF=/home/pato/Skrivebord/HPV16_projekt
SuperRunName=$SuperRunFolder
VirSupOut=$MainF/VirStrain_run/$SuperRunName
VirRunOut_run=$VirSupOut/$RunName

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder

###################### Create directory for each reference ##########################
# Finder navne for hver ref i References overmappe (eksl. undermapper)
find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_"${RunName}".txt
RefListOrigin=$(< $RefdF/RefSubtyper_"${RunName}".txt) # Finder navn på alle referencer i mappe

# HER DEFINERES REFERENCER
# RefListCalls=K02718.1_revised # Husk at comment RefList=$(< $RevRefs) og uncomment RefList=$RefListCalls
RefListCalls=$(< $VirRunOut_run/SubtypeCall.txt) # Uncomment hvis der skal bruges bedste call fra VirStrain som reference 



# Clearer RefList for referencer som ikke er i RefList mappe (eks fordi der ikke kunne kaldes op til 3 mulige typer og der derfor er blevet trukket forkert string ud)
rm -f $VirRunOut_run/Revised_SubTypeCalls.txt # Sletter for at sikre at der ikke appendes til fil fra tidligere test kørsel
touch $VirRunOut_run/Revised_SubTypeCalls.txt
RevRefs=$VirRunOut_run/Revised_SubTypeCalls.txt

echo Har kald: $RefListCalls

############## Clear unviable subtype calls #################
if [ "$customRef" = false ]; then
for refCall in $RefListCalls; do

	for refType in $RefListOrigin; do
	if [ "$refType" == "$refCall" ]; then
	Check=$(grep "$refCall" $VirRunOut_run/SubtypeCall.txt)
	if [ ${#Check} -gt 0 ]; then
	echo $refCall >> $RevRefs
	fi
	fi
	done

done
else 
	echo $RefListCalls >> $RevRefs
fi
##############################################################
