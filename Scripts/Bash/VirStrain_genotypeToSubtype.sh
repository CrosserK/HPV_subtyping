#!/bin/bash
set -e
set -u
set -o pipefail

# VirStrain subtyping

#!/bin/bash
set -e
set -u # Exit if undeclared variables

# VirStrain subtyping

# define help function
function help ()
{
	printf "Usage: %s: some var is missing 
	[-h|--help]\n" "$(basename "$0")" >&2
    exit
}

# Reading input parameters
while [[ "$#" -gt 0 ]]; do
	case "$1" in
        -r|--runname) RunName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-f|--fastq) FastQFile="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-s|--TopRunName) TopRunName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -m|--mainf) MainF="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -v|--virdb) VirStrain_db="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -c|--maincallfullname) MainCallFullName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

mkdir -p $MainF/VirStrain_run/$TopRunName/$RunName
mkdir -p $MainF/GenotypeCalls/$TopRunName/
#Define Folders and params
FQin=$MainF/FASTQ/${FastQFile}.fastq
SuperRunOut=$MainF/VirStrain_run/$TopRunName
Resultsout=$SuperRunOut/$RunName
GenoOut=$MainF/GenotypeCalls/$TopRunName/
VirStrain_db=$VirStrain_customdb

if [[ -d $VirStrain_db ]]; then

	cd ~/VirStrain

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
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $GenoOut/"${FastQFile}".txt
		grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $GenoOut/"${FastQFile}"_SplitTo.txt

	else
		echo $MainCallFullName > $Resultsout/SubtypeCall.txt
		echo "$MainCallFullName" > $GenoOut/${FastQFile}.txt
		echo "$MainCallFullName" > $GenoOut/${FastQFile}_SplitTo.txt
	fi

# Hvis ingen subtypedatabase, benyt da overtype	
else
	echo $MainCallFullName > $Resultsout/SubtypeCall.txt
	echo "(No subtype found)" >> $Resultsout/SubtypeCall.txt

	# Gemmer bedste strain til downstream scripts
	echo "$MainCallFullName" > $GenoOut/${FastQFile}.txt
	echo "$MainCallFullName" > $GenoOut/${FastQFile}_SplitTo.txt
	echo -e $FastQFile "$MainCallFullName" "(No Subtype found)" >> $SuperRunOut/VirStrain_summary.txt
fi
