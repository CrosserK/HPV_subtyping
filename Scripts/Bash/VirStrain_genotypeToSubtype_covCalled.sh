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
MainCallFullName=$6


##TEST
# RunName=pt_7.IonXpress_038_1
# FastQFile=pt_7.IonXpress_038 #input uden extenstion
# SuperRunName=KaroTest_1433_04102022
# MainF=/home/pato/Skrivebord/HPV_subtyping
# VirStrain_customdb=/home/pato/Skrivebord/HPV_subtyping/References/VirStrainDBs_HPV_sub/HPV16_VirStrainDB
# MainCallFullName=HPV16_K02718_1_revised
###

mkdir -p $MainF/VirStrain_run/$SuperRunName/"${RunName}"-"${MainCallFullName}"
#Define Folders and params
FQin=$MainF/FASTQ/${FastQFile}.fastq
SuperRunOut=$MainF/VirStrain_run/$SuperRunName
Resultsout=$SuperRunOut/"${RunName}"-"${MainCallFullName}"
VirStrain_db=$VirStrain_customdb

FQin=${FQin%.fastq}_filt.fastq

# Laver filer til downstream hvis de ikke eksisterer
if [[ ! -e $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt ]]; then
	touch $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
	touch $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo_sub.txt
fi


if [[ -d $VirStrain_db ]]; then

	cd ~/VirStrain

	echo Calling...
	echo python ~/VirStrain/VirStrain.py -i $FQin -d $VirStrain_db -o $Resultsout 

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

		# Checking if there is already a line with subtype, then concatenating new if true
		text=$(cat $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt)
		if [ ${#text} -gt 0 ]
		then
			old=$(awk 'FNR == 1 {print $1}' $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt)
			new=$(grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")')
			echo "$old"_"$new" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
		else
			grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
		fi

		# Also appending to end of line of _SplitTo file
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' >> $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo_sub.txt

	else
		echo $MainCallFullName > $Resultsout/SubtypeCall.txt
		# Checking if there is already a line with subtype, then concatenating new if true
		text=$(cat $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt)
		if [ ${#text} -gt 0 ]
		then
			old=$(awk 'FNR == 1 {print $1}' $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt)
			new=$(grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")')
			echo "$old"_"$new" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
		else
			grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
		fi		
		echo "$MainCallFullName" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt
		echo -e $FastQFile "$MainCallFullName" >> $SuperRunOut/VirStrain_summary.txt
	fi

# Hvis ingen subtypedatabase, benyt da overtype	
else
	echo $MainCallFullName > $Resultsout/SubtypeCall.txt

	# Gemmer bedste strain til downstream scripts
	text=$(cat $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt)
	if [ ${#text} -gt 0 ]
	then
		old=$(awk 'FNR == 1 {print $1}' $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt)
		echo "$old"_"$MainCallFullName" > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
	else
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt
	fi
	echo "$MainCallFullName" >> $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt
	echo -e $FastQFile "$MainCallFullName" >> $SuperRunOut/VirStrain_summary.txt
fi


mv $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_main.txt
mv $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo_main.txt

mv $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_sub.txt $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt
mv $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo_sub.txt $MainF/GenotypeCalls/$SuperRunName/${FastQFile}_SplitTo.txt


