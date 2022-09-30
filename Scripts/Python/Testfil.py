

# Vælger sekvens i range til at indsætte SNVs
subseq = "AAAAAAAAAAAAAAAAAAA"
length = len(subseq)

DNA = ""
for count in range(length):
    DNA += choice("AGTC")
subseq = DNA