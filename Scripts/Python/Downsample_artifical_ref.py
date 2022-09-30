import os

def open_ref_file(fasta):
    # Open fasta and return as fasta_string for global env
    fastafile = open(fasta, 'r')
    fasta_content = fastafile.read().replace("\n", "")
    fastafile.close()
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

def generate_deletion(seq, startpos, deletionlength):
    subseq = seq[startpos:startpos+deletionlength]
    flankbefore = seq[startpos - 10:startpos]
    flankafter = seq[startpos + deletionlength:startpos + deletionlength + 10]
    print("\nDeletion length: " + str(deletionlength) + "\n")
    newseq = str(seq[:startpos]) + str(
        seq[startpos + deletionlength:])  # Hvis insertion skal overwrite: +len(subseq) i sidste led
    return newseq

Choosefile = 1

if Choosefile == 1:
    name = "inversions_12_100.fa"
    amount = 1
elif Choosefile == 2:
    name = "chr17_insertions_12_100.fa"
    amount = 2
elif Choosefile == 3:
    name = "chr17_duplications_12_100.fa"
    amount = 2
else:
    print("Number not available")

oldfastafile = open_ref_file(name)
# Hvis duplications eller insertions er amount to delete, den del der er deleted fra deletions reference gange 2. Ved inversions er det gange 1
amountToDelete = sum([40, 100, 200, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000, 50000]) * 100 * amount

newseq = generate_deletion(oldfastafile, len(oldfastafile)-amountToDelete, amountToDelete)

if Choosefile == 1:
    newfastafile = open("inversions_d_12_100.fa", "w")
elif Choosefile == 2:
    newfastafile = open("chr17_insertions_d_12_100.fa", "w")
elif Choosefile == 3:
    newfastafile = open("chr17_duplications_d_12_100.fa", "w")
else:
    print("Number not available")


# Indsætter newlines efter hver 50. character korrekt igen
chrnumber=">chr17"
newseqheader=chrnumber
newseq = newseq.strip(chrnumber)
newseq = insert_newlines(newseq,50)
newseq = chrnumber + "\n" + newseq

# Writing to new reference fasta file
newfastafile.write(newseq)
newfastafile.close()