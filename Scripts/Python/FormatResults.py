import os, sys
import glob
import pandas as pd
from os import listdir
from os.path import isfile, join
from argparse import ArgumentParser

parser = ArgumentParser()

# define arguments. Argparse create help function from below help messages
parser.add_argument("-r", "--resultsF", required=True, help="Result folder for findhpvmutations script")
parser.add_argument("-t", "--topRunName", required=True, help="Result folder for findhpvmutations script")
parser.add_argument("-s","--sublineages", required=True, help='csv file of sublineages from PaVe. Format ex: Alpha-10,HPV11,A,A1,M14119')

# to allow passing lists, include the nargs='+'. The list must be given as ' ' seperated values, eg. --inlist a b c 
# to set a defualt, use default='', to set another variable, use dest='myargs' e.g.
# This will then be callable as args.myarg. To set a type, use type=int e.g.

args = parser.parse_args()

#Covarage matrix as csv file without qoutes, seperated by commas
resultsF = args.resultsF
topRunName = args.topRunName

indFile = "AnnotationIndividualFiles_"+topRunName+".txt"
sumFile = "AnnotationSummary_"+topRunName+".txt"

indFile = resultsF + "/" + indFile
sumFile = resultsF + "/" + sumFile

sublineagesFile = args.sublineages 


sumdf = pd.read_csv(sumFile, sep = "\t")
inddf = pd.read_csv(indFile, sep = "\t")
sublineages = pd.read_csv(sublineagesFile, sep = ",", header=None)

genes = list(sumdf['GENEID'].unique())
cols = ['Patient','Type','Reference matched','Sublineage']
cols.extend(genes)

# Removing error messages from table 
for item in cols.copy():
    if "DoesNotExist" in item:
        cols.remove(item)
for item in cols.copy():
    if "Vcf" in item:
        cols.remove(item)

# Changing genenames of ex gene-E4-172 to E4
newcols = []
for gene in cols:
    if "gene-" in gene:
        newgene = gene.split("-")[1]
        if newgene not in newcols:
            newcols.append(newgene)
    else:
        if gene not in newcols:
            newcols.append(gene)
cols = newcols

# Creating formatted table
sumJson = []
for patient in inddf.columns:
    # Making json with col names
    ptjson = {}
    for col in cols:
        ptjson.update({
            col : ""
        })
    
    # Getting patient name, reference and sublineage
    spl = patient.split(":")
    ref = spl[1]
    sublineage = ""
    type = ""
    # Getting sublineage class and type from the reference id
    
    subLineRow = 0
    for i in sublineages[4]:
        if i in ref:
            sublineage = sublineages[3][subLineRow]
            type = sublineages[1][subLineRow]
        subLineRow +=1

    ptjson.update({
        'Patient' : spl[0],
        'Type' : type,
        'Reference matched' : spl[1],
        'Sublineage' : sublineage
    })

    # For each row, get nuc and aa changes and match to gene key
    for row in inddf[patient]:
        # This check fails if value is NA
        if row == row:
            splt = row.split(" ")

            # If the cell contains a message such as "NoE4File", the column should be skipped. 
            if len(splt) == 2:
                continue

            if len(splt) == 4:

                nuc = splt[0]
                aa = splt[1]
                gene = splt[2]

                # Changing genenames of ex gene-E4-172 to E4
                if "gene-" in gene:
                    gene = gene.split("-")[1]
                
                if ptjson[gene] != "":
                    ptjson[gene].extend([[nuc,aa]])
                else:
                    ptjson.update({
                        gene : [[nuc,aa]]
                    })
            # Else if len is 5, there is a multi codon case. Because e.g c(36, 38) exists
            elif len(splt) == 5:
                nuc = splt[0]
                aa = str(splt[1:2]).replace(" ","")
                gene = splt[3]

                # Changing genenames of ex gene-E4-172 to E4
                if "gene-" in gene:
                    gene = gene.split("-")[1]
                
                if ptjson[gene] != "":
                    ptjson[gene].extend([[nuc,aa]])
                else:
                    ptjson.update({
                        gene : [[nuc,aa]]
                    })

    #jsonlist = [ptjson]
    sumJson.append(ptjson)

# Convert to dataframe and save to file
formDf = pd.DataFrame.from_dict(sumJson)
filename = resultsF+"/AnnotationFormatted_"+topRunName+".txt"
formDf.to_csv(filename, index = False)