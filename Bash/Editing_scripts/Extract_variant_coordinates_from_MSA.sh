	#!/bin/bash
set -e
set -u
set -o pipefail

# Finder SNP positioner fra multiple sequence fasta referencer og outputter som vcf. 

##### FREMGANG #####
# Lav multiple sequence alignment med MAFFT og output som fasta
# Sikrer at "main" strain er øverst i input (combined.fa)! SNP-sites benytter øverste som reference for resten af SNPs
# Husk at kopiere sublineage length ind i R script for header korrektion. Længder kommer ud nede i awk kommando
MainF=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
Input=HPV16_16Sub_rev.mafft
Output=HPV16_16Sub_relativeToK02718.1_revised.vcf
####################

snp-sites -v -o $MainF/$Output $MainF/$Input

# Get sublineage lengths, kopier den for "main strain" ind i R script for korrektion i header
# Dette skal gøres på unaligned file, da sekvensen ellers kan være forlænget med "-" pga gaps i alignment
Input1=${Input%mafft}fasta
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $MainF/$Input1


# VIGTIGT!!! Følg op med korrigering med Correct_MSA_vcf_File_to_main_strain.R script. 
