


# Fjern blanke linjer
sed '/^$/d' 

# Fjern alle RYSWKMBDHV fra nukleotid sekvenser
sed -e '/^>/! s/[RYSWKMBDHV]/N/g' /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_HG19_V5/HG19_HPV_V6.fasta /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_HG19_V5/HPV_HG19_V7.fasta 



sed -e '/^>/! s/[RYSWKMBDHV]/N/g' /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_HG19_V5/HG19_HPV_V6.fasta > temp.txt