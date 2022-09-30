
with open("chr17_duplications_12_100.fa", 'r') as f:
    content = f.read().replace("\n", "")
    print(content.index('gtggaatgaataaacaaatgacagcacattcgtgcaatggaatatgctgctgaggtgaaaatgaatgggtgcacacagcatgcattaacatgcataaatc')-6)

