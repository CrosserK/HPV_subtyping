#!/usr/bin/env python3

# Kør dette efter bash alignSubtypes til at generere gff3 fil:

import os
import sys
import pandas as pd
from argparse import ArgumentParser
# Åben gff3 fil

# arguments
parser = ArgumentParser()
parser.add_argument("-s", "--subtype", required=True,
                    help="hpvtype")
parser.add_argument("-m", "--maintype", required=True, help="output csv")

args = parser.parse_args()

subtype = args.subtype
maintype = args.maintype

#subtype = "KC470277_1"
#maintype = "HPV68"

gffFolder = "/home/pato/Skrivebord/HPV_subtyping/References/GFFfiles"
mainfolder = "/home/pato/Skrivebord/HPV_subtyping/Scripts/CreateDBs/E4_splice_update/Organised_Fastas/"+maintype
e4fixSavePath = "/home/pato/Skrivebord/HPV_subtyping/References/E4fixCoords"



# Checking that file does not already have AUHpred 
txt_file = open(gffFolder+"/"+subtype+".gff3", "r")
file_content = txt_file.read()
if "AUHpred" in file_content:
    sys.exit("File already have AUHpred")
else:
    print("updating file")

txt_file.close()


gfffile = gffFolder+"/"+subtype+".gff3"
#infile = pd.read_csv(gfffile)
with open(gfffile, "r") as file:
    gffstring = file.readlines()
    gffHeader1 = gffstring[0]
    gffHeader2 = gffstring[2]


# gfffile = gffFolder+subtype+".gff3"
# infile = pd.read_csv(gfffile)

# gffHeader1 = infile.columns[0]
# gffHeader2 = list(infile.loc[0])[0]
# # expanding file to create dataframe
# newlist = []

# for i in range(1,len(infile)):
#     newlist.append(list(infile.loc[i])[0])
# df = pd.DataFrame(newlist)
# fstSplit = df[0].str.split('\t', 8, expand=True)
# scnsplit = fstSplit[8].str.split(';', 8, expand=True)
# newdf = pd.concat([fstSplit,scnsplit], axis = 1)
# length = len(newdf.columns)
# for i in range(length):
#     newdf.columns.values[i] = i

coords = pd.read_csv(mainfolder+"/"+subtype+"_E4coords.csv", header = None, sep = " ")

# Tæller længde fra cigar string
from itertools import groupby

def query_len(cigar_string):
    """
    Given a CIGAR string, return the number of bases included in the
    query sequence.
    """
    read_consuming_ops = ("M", "I", "S", "=", "X")
    read_subtracting_ops = "D"
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            result += length
        if op in read_subtracting_ops:
            result -= length
    return result

e41len = query_len(str(coords.iloc[0,1]))
e42len = query_len(str(coords.iloc[1,1]))

e4s_1 = coords.iloc[0,0]
e4e_1 = e4s_1+e41len-1
e4s_2 = coords.iloc[1,0]
e4e_2 = e4s_2+e42len-1

# update the file
fs = gffHeader2.split("\t")
# Get the length of the first splice section to infer CDS phase. +1 because the positions are inclusive.
cdsPhaseCalc = (e4s_2 - e4s_1 + 1) % 3

if cdsPhaseCalc == 0:
    cdsPhase = 0
elif cdsPhaseCalc == 1:
    cdsPhase = 2
elif cdsPhaseCalc == 2:
    cdsPhase = 1
else:
    sys.exit("CDSphase does not go into 0,1 or 2")

row1 = [fs[0],"AUHpred","region",fs[3],fs[4],fs[5],fs[6],fs[7],fs[8]]
row2 = [fs[0],"AUHpred","gene",e4s_1,e4e_2,'.','+','.',"ID=E4_splice"]
row3 = [fs[0],"AUHpred","mRNA",e4s_1,e4e_2,'.','+','.',"ID=E4_spliceRNA;Parent=E4_splice"]
row4 = [fs[0],"AUHpred","CDS",e4s_1,e4e_1,'.','+','0',"ID=E4cds1;Parent=E4_spliceRNA"]
row5 = [fs[0],"AUHpred","CDS",e4s_2,e4e_2,'.','+',cdsPhase,"ID=E4cds2;Parent=E4_spliceRNA"]

list_lists=[row1,row2,row3,row4,row5]
spliceDf = pd.DataFrame(list_lists)
spliceDf.iloc[0,8] = spliceDf.iloc[0,8].split(":")[0]

# save to correct position
# save to correct position
# Add header and save to file

row0 = "##gff-version 3"
row01 = gffHeader1

file1 = open(e4fixSavePath+"/"+subtype+".gff3", "w")
file1.writelines(row0+"\n")
file1.writelines(row01+"\n")
file1.close()

spliceDf.to_csv(e4fixSavePath+"/"+subtype+".gff3", header = None, index = None, sep = "\t", mode = "a")
# Now update the original gff3 file

# Remove existing E4 lines

# Add new ones
row1 = [fs[0],"AUHpred","gene",e4s_1,e4e_2,'.','+','.',"ID=E4_splice"]
row2 = [fs[0],"AUHpred","mRNA",e4s_1,e4e_1,'.','+','.',"ID=E4_spliceRNA1;Parent=E4_splice"]
row3 = [fs[0],"AUHpred","mRNA",e4s_2,e4e_2,'.','+','.',"ID=E4_spliceRNA2;Parent=E4_splice"]
row4 = [fs[0],"AUHpred","CDS",e4s_1,e4e_1,'.','+','0',"ID=E4exon1;Parent=E4_spliceRNA1"]
row5 = [fs[0],"AUHpred","CDS",e4s_2,e4e_2,'.','+',cdsPhase,"ID=E4exon2;Parent=E4_spliceRNA2"]

list_lists=[row1,row2,row3,row4,row5]
spliceDf = pd.DataFrame(list_lists)
spliceDf.iloc[0,8] = spliceDf.iloc[0,8].split(":")[0]

spliceDf.to_csv(gffFolder+"/"+subtype+".gff3", header = None, index = None, sep = "\t", mode = "a")
