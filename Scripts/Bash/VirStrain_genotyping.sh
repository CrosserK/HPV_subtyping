#!/bin/bash
set -e
set -u

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
		-f|--fastqaddr) FQAddr="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-s|--TopRunName) TopRunName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -m|--mainf) MainF="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -v|--virdb) VirStrain_db="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -t|--typetier) typeTier="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -g|--genotypetop) genoTypeTop="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

FQName=$(basename $FQAddr .fastq)

mkdir -p $MainF/VirStrain_run/$TopRunName/GenotypeCalls/$FQName
mkdir -p $MainF/GenotypeCalls/$TopRunName/${FQName}/
#Define Folders and params
FQin=$MainF/FASTQ/${FQName}.fastq
SuperRunOut=$MainF/VirStrain_run/$TopRunName/GenotypeCalls
Resultsout=$SuperRunOut/${FQName}

cd ~/VirStrain

# python VirStrain_build.py -i '/home/pato/Skrivebord/HPV16_projekt/References_andre/16substrain_HPV16_sublineages.mafft' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_16_virstrain' -s 1
python /home/pato/VirStrain/VirStrain.py -i $FQin -d $VirStrain_db -o $Resultsout 

cd

# Finder kaldte subtyper:
# Printer første felt af 2,3 og 4. linje og fjerner ">" (1. linje er header i rapport)
grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $Resultsout/SubtypeCall.txt
grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 3 {print $1}' | awk 'sub(/^>/, "")' >> $Resultsout/SubtypeCall.txt
grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 4 {print $1}' | awk 'sub(/^>/, "")' >> $Resultsout/SubtypeCall.txt

echo -e $FQName '\t' "$(grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")')" >> $SuperRunOut/VirStrain_summary.txt
echo -e "$(grep -A10 Top10_Score_Strains ${Resultsout}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")')" >> $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt