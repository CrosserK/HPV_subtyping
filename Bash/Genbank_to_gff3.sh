


bp_genbank2gff 'in.gb' --viral --stdout > 'out.gff3'
# Make sure the name is correct. The script seems to remove the ".1" from "K02718.1" in the gff3 file

# HUSK følgende 3 korrigeringer!
# Lav korrekt navn, eks:
sed 's/K02718/K02718.1/g'

# Fjern .t00:
sed 's/.t00//g'

# Fjern _ for at R script ikke kommer til at separere forkert:
sed 's/_ALPHA/-ALPHA/g'

# Sæt CDS phases til 0 (Tjek om de er 0 i genbank fil), byg videre på:
# grep "CDS" $file | awk "CDS" {print $8} 


