
# Search for subsequences in fasta file

file=

seqkit locate $file -p "CCCCATGCACCAAAAAAACTAAG" -m 3 
# -m defines number of allowed mismatches