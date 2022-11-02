import os
from random import choice

# Open fasta reference and return it as string
def open_ref_file(fasta):
    # Open fasta and return as fasta_string for global env
    fastafile = open(fasta, 'r')
    fasta_content = fastafile.read().replace("\n", "")
    return str(fasta_content)

def insert_newlines(string, every=50):
    # Indsætter newline efter hver x'ne character hvor "every" er x
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

# Insert chr number. It will return the position of end of the written string
def find_chr_pos(fasta, chr):
    start_index = fasta.find(chr) + 1
    end_index = start_index + len(chr) - 1
    return end_index

# Inverts the character order of a chosen sequence from start to end position and
# returns newseq with full reference and inversion
def generate_SNVs_in_range(seq, startpos, length,SNVpercent):
    endpos = startpos+length

    # Vælger sekvens i range til at indsætte SNVs
    subseq = seq[startpos:endpos]

    # This function generates an inversion of the input string
    # print("\nFull sequence length: " + str(len(subseq)) + "\n")
    # Printing slice sequence
    # print("Chosen subsequence:\n" + str(subseq))
    # String splitter
    # Inversion of slice
    # print("\nInverted subsequence:\n" + invertedsubseq)

    SV_seq_details_file.write("\nStartpos: " + str(startpos-6) + "\n" + "Full sequence length: " + str(len(subseq)) + "\n" + "----Chosen subsequence----\n" + str(subseq) + "\n" + "----Inverted subsequence----\n" + invertedsubseq)
    newseq = str(seq[:startpos]) + str(invertedsubseq) + str(seq[endpos:])

    # Sørger for at indsætte ny sekvens som overskriver gammel
    # newseq = str(seq[:startpos]) + str(subseq) + str(seq[startpos+length:])
    return newseq

# Running the program:
# Open fasta as fastafile
fastafile = open_ref_file("chr17_control.fa")  # Opens ref as string

SVlengths = [40, 100, 200, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000, 50000] * 100

# Laver koordinater til SV med 75000 mellemrum
startpos=500000
newlist=[startpos]
# Definerer variabel som integer:
pos=startpos
stdspacer = 25000

for tenrepeats in range(0, len(SVlengths)): # laver 120 fra 0 til 119
    newrepeatpos=pos
    spacer = stdspacer + SVlengths[tenrepeats]
    pos = newrepeatpos + spacer
    if pos > 21782000:
        pos = pos + 5155000
    print(pos)
    newlist.append(pos)
    # resetter lige position til udgangspunkt
    if pos > 21782000:  #Dodger NNN region
        pos = pos - 5155000
# Fjerner sidste punkt som blev sat på fordi loop laver 120 tal, men startpos allerede er sat ind (så der kommer 121)
newlist.pop()
startcoordinates=newlist

headerline = ">chr17"
newseq = fastafile
SV_infofile= open("SV_info_file_chr17_inversions_12_100.txt","a") # a option for append
SV_infofile.write("SVs inserted at: \n")

for x in startcoordinates:
    y = startcoordinates.index(x)
    InversionLength = SVlengths[y]

    # Insert chr number. It will return the position of end of the written string
    chrstart = find_chr_pos(fastafile, headerline)

    # korrigerer position
    SVposition = chrstart + x

    # Generating SV's, inversion stores as "newseq" globally
    SV_seq_details_file = open("SV_seq_details_chr17_inversions_12_100.txt", "a")  # a option for append
    newseq = generate_inversion(newseq, SVposition, InversionLength)
    SV_seq_details_file.close()
    infotext="\nStartcoordinate: " + str(x) + " " + "Length: " + str(InversionLength)
    print(infotext)
    SV_infofile.write(infotext)

SV_infofile.close()

# Formatterer infofil
with open("SV_seq_details_chr17_inversions_12_100.txt") as infile, open("SV_seq_details_chr17_inversions_12_100_formatted.txt", 'w') as outfile:
  for line in infile:
    if len(line) > 50:
      outfile.write('\n'.join(line[i:i+50] for i in range(0,len(line), 50)))
    else:
      outfile.write(line)
os.remove("SV_seq_details_chr17_inversions_12_100.txt")

# Creating new fasta file with SV's
newfastafile = open("chr17_inversions_12_100.fa", "w")

# Indsætter newlines efter hver 50. character korrekt igen
newseqheader=headerline
newseq = newseq.strip(headerline)
newseq = insert_newlines(newseq,50)
newseq = headerline + "\n" + newseq

# Writing to new reference fasta file
newfastafile.write(newseq)
newfastafile.close()

# Hvis SV'er opstår i NN regioner (søg i seq detail fil) indsæt da nye conditions i NN dodger loop

