#!/bin/bash

# HPV typer er splittet op efter nummer med alle undertype i mappe herunder.
# Der skal nu alignes til fundne sekvenser 1 overtype af gangen. De fundne sekvenser findes fra createE4Fixsubtype.ipynb, 
# ved at indtaste hovedtype E4 sekvenser fra PaVE. Der er 28 hvor det skal gøres med. 
# Nedenstående er udvikling af automatisk alignment og finde ud af om de giver mening, for hver HPV gruppe. 


# Reading input parameters
while [[ "$#" -gt 0 ]]; do
	case "$1" in
        -i|--input) hpvName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
	   -f|--location) location="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

#hpvName=HPV11

location=/home/pato/Skrivebord/HPV_subtyping/Scripts/CreateDBs/E4_splice_update/Organised_Fastas/$hpvName

#cd $location

fastq=${location}/${hpvName}_E4_seq.fastq


# Keep up on errors:
log=()





# Indexerer hver fasta fil som ikke er maintype
for subtype in ${location}/*.fasta; do
	BN=$(basename $subtype)

	# Testing that the basename does not contain the hpvname to avoid getting the main hpv type
	if [[ $BN != *"${hpvName}"* ]]; then

		# Kører nu algoritmer for at sammenligne E4 sekvenser fra mainline til hver subtype som ikke er mainline

		echo Processing subtype: "$(basename "${subtype%.fasta}")"
		bwa index $subtype
		fastq=${location}/${hpvName}_E4_seq.fastq
		#bwa aln -i 1 -m 2 -o 4 $subtype $fastq > ${subtype%.fasta}.sam
		bwa mem -T 9 -k 3 -L 15 $subtype $fastq > ${subtype%.fasta}.sam
		samtools view -Sb ${subtype%.fasta}.sam > ${subtype%.fasta}.bam; 
		samtools sort ${subtype%.fasta}.sam -o ${subtype%.fasta}.sorted.bam
		samtools index ${subtype%.fasta}.sorted.bam

		#Check that they start with start codon and end with stop codon
		file=${subtype%.fasta}.sorted.bam
		# Getting first codon of 1st splice
		fstCodon=$(samtools view $file | awk '{ print $10 }' | head -n 1 | cut -c 1-3)
		# Getting last codon of 2nd splice
		lstCodon=$(samtools view $file | awk '{ print $10 }' | head -n 2 | tail -n 1 | rev | cut -c 1-3 | rev)
		
		# Check that 1st codon is start codon
		if [ $fstCodon = ATG ]; then
		cod1=pass
		else
		cod1=fail
		# Add new element at the end of the array
		arrVar+=(${subtype})
		arrVar+=("Did not have start codon")
		#Output to errors file
		fi

		# Check that last codon is stop codon
		if [ $lstCodon = TAA ] || [ $lstCodon = TAG ] || [ $lstCodon = TGA ]; then
		cod2=pass
		else
		cod2=fail
		arrVar+=(${subtype})
		arrVar+=("Did not have stop codon")
		#Output to errors file
		fi

		if [ $cod1 = pass ] && [ $cod2 = pass ]; then
		echo saving "$(basename "${subtype%.fasta}")" to file
		# If passed, save to file 
		samtools view $file | awk '{ print $4 " " $6 }' > "${subtype%.fasta}"_E4coords.csv
		
		# creating E4 splice gff3 files
		python createGffForSubtype.py -s "$(basename "${subtype%.fasta}")" -m ${hpvName}

		else
			echo File for "$(basename "${subtype%.fasta}")" not created due to above error
			sleep 5
		fi


		# Kør python script til dette

	
	fi

	sleep 2

done

# Iterate the loop to read and print each array element
for value in "${log[@]}"
do
     echo $value
done


