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
		-f|--fastqaddr) FQAddr="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-s|--TopRunName) TopRunName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -m|--mainf) MainF="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -v|--virdb) VirStrain_db="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -c|--maincallfullname) MainCallFullName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -t|--typetier) typeTier="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
		-g|--subtypetop) subTypeTop="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

# Incrementing input typeTier to output typeTier
let "typeTier=typeTier+1"

#Define Folders and params
FQName=$(basename $FQAddr .fastq)

VirTopFOut=$MainF/VirStrain_run/$TopRunName
VirRunOut=$VirTopFOut/"${FQName}"-"${MainCallFullName}"
mkdir -p $VirRunOut

subOut=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}.txt 
subOutSplit=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt


# Making downstream files if it is the first time this script is run
if [[ ! -e $subOut ]]; then
	touch $subOut
	touch $subOutSplit
fi

# Finding subtype
if [[ -d $VirStrain_db ]]; then
	cd ~/VirStrain # VirStrain require you to be in the folder to run. 
	python /home/pato/VirStrain/VirStrain.py -i ${FQAddr} -d $VirStrain_db -o $VirRunOut 
	cd

	# Check if a subtype has been found
	if [ -f ${VirRunOut}/VirStrain_report.txt ]; 
	then	
		# Print top 3. (1. line is header in report)
		grep -A10 Top10_Score_Strains ${VirRunOut}/VirStrain_report.txt | awk 'FNR == 2 {print $1}' | awk 'sub(/^>/, "")' > $VirRunOut/SubtypeCall.txt
		grep -A10 Top10_Score_Strains ${VirRunOut}/VirStrain_report.txt | awk 'FNR == 3 {print $1}' | awk 'sub(/^>/, "")' >> $VirRunOut/SubtypeCall.txt
		grep -A10 Top10_Score_Strains ${VirRunOut}/VirStrain_report.txt | awk 'FNR == 4 {print $1}' | awk 'sub(/^>/, "")' >> $VirRunOut/SubtypeCall.txt
		# Print top 3 to collective file for run
		top1=$(awk 'FNR == 1 {print $1}' $VirRunOut/SubtypeCall.txt) 
		top2=$(awk 'FNR == 2 {print $1}' $VirRunOut/SubtypeCall.txt) 
		top3=$(awk 'FNR == 3 {print $1}' $VirRunOut/SubtypeCall.txt)
		echo -e ${FQName} '\t' $top1 '\t' $top2 '\t' $top3 >> $VirTopFOut/VirStrain_summary.txt

		# Now sending to GenotypeCalls folder, for downstream processing
		# Checking if there is already a line with subtype (this script has already been run for another genotype), then concatenating new one if true
		# Now will only look at highest scoring subtype for each maintype
		text=$(cat $subOut)
		if [ ${#text} -gt 0 ]
		then
			old=$(awk 'FNR == 1 {print $1}' $subOut)
			echo "$old"_"$top1" > $subOut
		else
			echo $top1 > $subOut
		fi

		# Also appending to end of line of _SplitTo file
		echo $top1 >> $subOutSplit

	else
		# Else if no subtype found, set subtype as genotype
		echo $MainCallFullName > $VirRunOut/SubtypeCall.txt
		top1=$(awk 'FNR == 1 {print $1}' $VirRunOut/SubtypeCall.txt) 
		# Checking if there is already a line with subtype (this script has already been run for another genotype), then concatenating new one if true
		# Now will only look at highest scoring subtype for each maintype
		text=$(cat $subOut)
		if [ ${#text} -gt 0 ]
		then
			old=$(awk 'FNR == 1 {print $1}' $subOut)
			echo "$old"_"$top1" > $subOut
		else
			echo $top1 > $subOut
		fi

		# Also appending to end of line of _splitTo file
		echo $top1 >> $subOutSplit

		echo -e $FQName "$MainCallFullName" >> $VirTopFOut/VirStrain_summary.txt
	fi

# If no subtypedatabase, use genotype as subtype. Subtype files are used for further processing. 
else
	echo $MainCallFullName > $VirRunOut/SubtypeCall.txt
	# Saving best strain
	# Checking if script has already been run, so there is already data in this text file
	text=$(cat $subOut)
	if [ ${#text} -gt 0 ]
	then
		old=$(awk 'FNR == 1 {print $1}' $subOut)
		echo "$old"_"$top1" > $subOut
	else
		echo $top1 > $subOut
	fi
	echo "$MainCallFullName" >> $subOutSplit
	echo -e $FQName "$MainCallFullName" >> $VirTopFOut/VirStrain_summary.txt
fi


