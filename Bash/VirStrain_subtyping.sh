#!/bin/bash
set -e
set -u
set -o pipefail

# VirStrain subtyping

RunName=$1
FastQFile=$2 #input uden extenstion
VirStrain_customdb=$3
cutadaptMinsize=$4
cutadaptMaxsize=$5

QualTrim=20
#RunName=Pt_11_DNA.IonXpress_085_run
#FastQFile=Pt_11_DNA.IonXpress_085 #input uden extenstion
#VirStrain_customdb=HPV16_16_virstrain

###################### EDIT ##########################
#Define Folders and params
MainF=/home/pato/Skrivebord/HPV16_projekt
FQin=$MainF/FASTQ/${FastQFile}.fastq
Resultsout=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/$RunName
VirStrain_db=$MainF/References/$VirStrain_customdb

# Kvalitetesikrer
cutadapt -q $QualTrim -m $cutadaptMinsize -M $cutadaptMaxsize $FQin -o ${FQin%.fastq}_filt.fastq

FQin=${FQin%.fastq}_filt.fastq

cd VirStrain

# python VirStrain_build.py -i '/home/pato/Skrivebord/HPV16_projekt/References_andre/16substrain_HPV16_sublineages.mafft' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_16_virstrain' -s 1
python /home/pato/VirStrain/VirStrain.py -i $FQin -d $VirStrain_db -o $Resultsout 

cd

# Samler resultat
echo -e $FastQFile '\t' $(cat ${Resultsout}/VirStrain_report.txt | sed 's/^.*;QD=\(>*\);.*$/\1/' | awk 'NR==3' | cut -f3 | sed 's/>//' ) >> $MainF/VirStrain_run/VirStrain_summary.txt