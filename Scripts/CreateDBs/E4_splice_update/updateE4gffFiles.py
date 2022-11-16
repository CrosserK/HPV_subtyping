import os
import sys
from Bio import SeqIO
import glob
import subprocess

# Insert
maintype = "HPV53"
e4_1s = 892
e4_1e = 907
e4_2s = 3343
e4_2e = 3653
mainfolder = "Organised_Fastas/"+maintype

# find maintype fasta
fasta = glob.glob(mainfolder+"/"+maintype+"*.fasta")
# Check that there is only one match:
if len(fasta) != 1:
    sys.exit("STOP, more than one match for mainline fasta")
else:
    fasta = fasta[0]
item = SeqIO.parse(fasta, "fasta")
# getting record
for i in item:
    record = i

e4_1 = record.seq[e4_1s-1:e4_1e]
e4_2 = record.seq[e4_2s-1:e4_2e]


seqLen1 = len(e4_1)
seqLen2 = len(e4_2)

with open(mainfolder+"/"+maintype+"_E4_seq.fastq", "w") as text_file:
    text_file.write("@E4_splice1\n")
    text_file.write(str(e4_1)+"\n")
    text_file.write("+\n")
    text_file.write("~"*seqLen1+"\n")
    
    text_file.write("@E4_splice2\n")
    text_file.write(str(e4_2)+"\n")
    text_file.write("+\n")
    text_file.write("~"*seqLen2+"\n")

# Kør bash script for aligning og codon check
val = subprocess.check_call("./alignSubtypes.sh -i '%s' -f '%s'" % (maintype, mainfolder), shell=True)