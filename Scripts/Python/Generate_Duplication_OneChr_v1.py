
import os

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

# Inserts the character order of a chosen sequence from start position and
# returns newseq with full reference and insertions of different lengths of the given sequence (repeating it if longer than given)
def generate_duplication(seq, startposOfDup, length):
    lenToDup = length # Sørger for at hele lister bliver brudt op i enkelt værdier inden de sendes videre i programmet
    seqToDup = str(seq[startposOfDup: startposOfDup+lenToDup])
    print("\nFull duplication length: " + str(len(seqToDup)) + "\n")
    # Printing slice sequence
    print("Chosen duplication:\n" + str(seqToDup))
    flankbefore = seq[startposOfDup - 10:startposOfDup]
    flankafter = seq[startposOfDup+lenToDup:startposOfDup+lenToDup + 10]
    print("Duplication length: " + str(len(seqToDup)) + "\n")
    SV_seq_details_file.write("\nStartpos: " + str(startposOfDup - 6) + "\n" + "Full sequence length: " + str(
        len(seqToDup)) + "\n" + "----Chosen subsequence duplicated----\n" + str(
        seqToDup) + "\n" + "-----Flanked by-----\n" + flankbefore + " & " + flankafter + "\n")
    newseq = str(seq[:startposOfDup]) + str(seqToDup) + str(seq[startposOfDup:]) # Hvis insertion skal overwrite: +len(seqToDup) i sidste led
    return newseq


# Running the program:
# Open fasta as fastafile
fastafile = open_ref_file("chr17_control.fa")  # Opens ref as string

# chr17= 84922596bp. Mellem 70000 og 83.247.450 udelader start og slut N'er (10.118.250)
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
    # resetter position til udgangspunkt
    if pos > 21782000:
        pos = pos - 5155000

# Fjerner sidste punkt som blev sat på fordi loop laver 120 tal, men startpos allerede er sat ind (så der kommer 121)
newlist.pop()
startcoordinates=newlist

#for tenrepeats in range(0,10): # laver 10 fra 0 til 9
#    newrepeatpos=pos
#    for numberofSV in range (1,13): # Laver 12, fra 1 til 12.
#        pos = newrepeatpos + spacer * numberofSV
#        print(pos)
#        newlist.append(pos)
# Fjerner sidste punkt som blev sat på fordi loop laver 120 tal, men startpos allerede er sat ind (så der kommer 121)
#newlist.pop()
#startcoordinates=newlist


chrnumber = ">chr17"
newseq = fastafile
SV_infofile= open("SV_info_file_chr17_duplications_12_100.txt","a") # a option for append
SV_infofile.write("SVs inserted at: \n")

for x in startcoordinates:
    y = startcoordinates.index(x)
    InsertionLength = SVlengths[y]

    # Insert chr number. It will return the position of end of the written string
    chrstart = find_chr_pos(fastafile, chrnumber)
    # korrigerer position pga sekvense ændres løbende af de deleterede områder
    if y > 0:
        SVposition = chrstart + x + sum(SVlengths[:y])
    else:
        SVposition = chrstart + x


    # Generating SV's, inversion stores as "newseq" globally
    SV_seq_details_file = open("SV_seq_details_chr17_duplications_12_100.txt", "a")  # a option for append
    newseq = generate_duplication(newseq, SVposition, InsertionLength)
    SV_seq_details_file.close()
    infotext="\nStartcoordinate: " + str(x) + " " + "Length: " + str(InsertionLength)
    print(infotext)
    SV_infofile.write(infotext)

SV_infofile.close()

# Formatterer infofil
with open("SV_seq_details_chr17_duplications_12_100.txt") as infile, open("SV_seq_details_chr17_duplications_12_100_formatted.txt", 'w') as outfile:
  for line in infile:
    if len(line) > 50:
      outfile.write('\n'.join(line[i:i+50] for i in range(0,len(line), 50)))
    else:
      outfile.write(line)
os.remove("SV_seq_details_chr17_duplications_12_100.txt")

# Creating new fasta file with SV's
newfastafile = open("chr17_duplications_12_100.fa", "w")

# Indsætter newlines efter hver 50. character korrekt igen
newseqheader=chrnumber
newseq = newseq.strip(chrnumber)
newseq = insert_newlines(newseq,50)
newseq = chrnumber + "\n" + newseq

# Writing to new reference fasta file
newfastafile.write(newseq)
newfastafile.close()

# Hvis SV'er opstår i NN regioner (søg i seq detail fil) indsæt da nye conditions i NN dodger loop

