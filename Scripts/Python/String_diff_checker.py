# Checks if two text files differ and prints the differing line

f1=open("chr17.fa","r")
f2=open("chr17_instertedSVs.fa","r")
for line1 in f1:
    for line2 in f2:
        if line1 != line2:
            print("File1:\n" + line1 + "File2:\n" + line2)
        break
f1.close()
f2.close()