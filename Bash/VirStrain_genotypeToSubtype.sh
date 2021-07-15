#!/bin/bash
set -e
set -u
set -o pipefail

# VirStrain subtyping

RunName=$1
FastQFile=$2 #input uden extenstion
SuperRunName=$3
MainF=$4
VirStrain_customdb=$5
FoundMainType=$6
MainCallFullName=$7

##TEST
#RunName=Pt_33_RNA.IonXpress_087_run
#FastQFile=Pt_33_RNA.IonXpress_087 #input uden extenstion
#SuperRunName=HPV_double_subtyping_1410_29062021
#MainF=/home/pato/Skrivebord/HPV16_projekt
#VirStrain_customdb=/home/pato/Skrivebord/HPV16_projekt/References/VirStrainDBs_HPV_sub/HPV16_VirStrainDB
#FoundMainType=HPV16
###

mkdir -p $MainF/VirStrain_run/$SuperRunName/$RunName
#Define Folders and params
FQin=$MainF/FASTQ/${FastQFile}.fastq
SuperRunOut=$MainF/VirStrain_run/$SuperRunName
Resultsout=$SuperRunOut/$RunName
VirStrain_db=$VirStrain_customdb

FQin=${FQin%.fastq}_filt.fastq

if [[ -d $VirStrain_db ]]; then

	cd VirStrain

	# python VirStrain_build.py -i '/home/pato/Skrivebord/HPV16_projekt/References_andre/16substrain_HPV16_sublineages.mafft' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_16_virstrain' -s 1
	python /home/pato/VirStrain/VirStrain.py -i $FQin -d $VirStrain_db -o $Resultsout 

	cd

	# Checker om der blev fundet en subtype og om der derfor blev lavet en rapport fil
	if [ -f ${Resultsout}/VirStrain_report.txt ]; then
		
		# Finder kaldte subtyper:
		# Printer første felt af 2,3 og 4. linje og fjerner ">" (1. linje er header i rapport)
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $Resultsout/SubtypeCall.txt
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 3 {print $1}' | awk 'sub(/^>/, "")' >> $Resultsout/SubtypeCall.txt
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 4 {print $1}' | awk 'sub(/^>/, "")' >> $Resultsout/SubtypeCall.txt

		echo -e $FastQFile '\t' $(grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")') >> $SuperRunOut/VirStrain_summary.txt

		# Gemmer bedste strain til downstream scripts
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt

	else
		echo $MainCallFullName > $Resultsout/SubtypeCall.txt
		echo "(No subtype references)" >> $Resultsout/SubtypeCall.txt
		echo "$MainCallFullName" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt
		echo "$MainCallFullName" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt
		echo -e $FastQFile "$MainCallFullName" "(No Subtype found)" >> $SuperRunOut/VirStrain_summary.txt
	fi

# Hvis ingen subtypedatabase, benyt da overtype	
else
	echo $MainCallFullName > $Resultsout/SubtypeCall.txt
	echo "(No subtype found)" >> $Resultsout/SubtypeCall.txt

	# Gemmer bedste strain til downstream scripts
	echo "$MainCallFullName" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt
	echo "$MainCallFullName" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt
	echo -e $FastQFile "$MainCallFullName" "(No Subtype found)" >> $SuperRunOut/VirStrain_summary.txt
fi
