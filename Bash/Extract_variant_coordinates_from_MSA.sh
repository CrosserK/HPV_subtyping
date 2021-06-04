#!/bin/bash
set -e
set -u
set -o pipefail

# Finder SNP positioner fra multiple sequence fasta referencer og output som vcf. 

##### FREMGANG #####
# Lav multiple sequence alignment med Clustal Omega Online tool eller MAFFT og output som fasta
# Sikrer at "main" strain er øverst i input (combined.fa)! SNP-sites benytter øverste til at kigge på resten af SNPs med
# Husk at kopiere sublineage length ind i R script for header korrektion. Længder kommer nede i awk
MainF=/home/pato/Skrivebord/HPV16_projekt/References_andre
# Kombinerer alle fasta i folder

Input=15_HPV_sublineages.mafft
Output=15_HPV_sublineages_aligned_MAFFT_relativeToK02718.1.vcf
####################


# Sammensætte fasta'er til Clustal Omega
#cd $MainF/References
#cat *.fasta > $MainF/combined.fasta
#cd

snp-sites -v -o $MainF/$Output $MainF/$Input

# Get sublineage lengths, kopier den for "main strain" ind i R script for korrektion i header
# Dette skal gøres på unaligned file, da sekvensen ellers kan være forlænget med "-" pga gaps i alignment
Input=15_HPV_sublineages.fasta
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $MainF/$Input

# Følg op med korrigering med Correct_MSA_vcf_File_to_main_strain.R


######## GAMMELT ###############


TargetFolder='/home/pato/Skrivebord/HPV16_projekt'

cd $TargetFolder/ReferenceDetails

awk 'NR % 12 == 0' Clustal_HPV_all.clustal > stars.txt
cat stars.txt | sed 's/^.\{16\}//g'| sed ':a; N; $!ba; s/\n//g' > stars1.txt # Removes the first 16 characters (incl whitespace), then removes all spaces
cat stars1.txt | grep -aob ' ' | grep -oE '[0-9]+' |  awk '{print $1+1}' > coordinates.txt #Finds position of whitespaces and adds 1 to give exact position

# OBS ikke færdig, den tæller lidt for langt!


# Multiple sequence alignment: (input ren fasta fil lavet med kalign combine men uden headeren)
kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa

kalign '/home/pato/Skrivebord/HPV16_projekt/References/K02718.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/AF536179.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/AF536180.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/AF472509.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/AF402678.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/AY686579.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/HQ644298.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/AF534061.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/HQ644236.1.fasta' '/home/pato/Skrivebord/HPV16_projekt/References/HQ644257.1.fasta' -f fasta > '/home/pato/Skrivebord/HPV16_projekt/References/HPV_all_aligned.afa'



# Find number of occurences
var="ACAGTCGATCG--CAGCTAG-ACA"
res="${var//[^-]}"
echo "${#res}"

#Find position of - in reference
