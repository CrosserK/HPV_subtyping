#!/bin/bash
set -e
set -u
set -o pipefail

##### FREMGANG #####
# Lav multiple sequence alignment med Clustal Omega Online tool og output som fasta
# Uncomment section hvis der skal sammensættes mange fasta
# Sikre at "main" strain øverst i input (combined.fa)! SNP-sites benytter øverste til at kigge på resten af SNPs med
# Find SNP positioner fra combined fasta referencer og output som vcf. 
# Husk at kopiere sublineage length ind i R script for header korrektion. Længder kommer nede i awk
MainF=/home/pato/Skrivebord/HPV16_projekt
# Kombinerer alle fasta i folder

Input=15_HPV_sublineages_aligned.afa
Output=15_HPV_sublineages_aligned.vcf
####################


# Sammensætte fasta'er til Clustal Omega
#cd $MainF/References
#cat *.fasta > $MainF/combined.fasta
#cd

# Get sublineage lengths, kopier dette ind i R script for korrektion i header
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $MainF/$Input

snp-sites -v -o $MainF/$Output $MainF/$Input






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
