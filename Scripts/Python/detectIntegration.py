import os, sys
import glob
import pandas as pd
from os import listdir
from os.path import isfile, join
from argparse import ArgumentParser

parser = ArgumentParser()

# define arguments. Argparse create help function from below help messages
# define arguments. Argparse create help function from below help messages
parser.add_argument("-i", "--infile", required=True, help="Coverage matrix as csv file without qoutes, seperated by commas")
parser.add_argument("-r", "--ref", required=True, help="Location of all valid references that can be linked to")
parser.add_argument("-o","--outpath", required=True, help='Location of where to save output files')
parser.add_argument("-g","--houseGeneEndRow", required=True, help='Last row containing housegenes. They should be from row 1 to this row')
parser.add_argument("-a","--amplPos", required=True, help='bed file with amplicon positions. See manual for format')

# to allow passing lists, include the nargs='+'. The list must be given as ' ' seperated values, eg. --inlist a b c 
# to set a defualt, use default='', to set another variable, use dest='myargs' e.g.
# This will then be callable as args.myarg. To set a type, use type=int e.g.

args = parser.parse_args()

amplPos = args.amplPos
houseGeneEndRow = args.houseGeneEndRow
#Covarage matrix as csv file without quotes, seperated by commas
file = args.infile

# Location of all valid reference names that can be given
path = args.ref

# Location of where to save output files
refSavePath = args.outpath


faF = [f for f in listdir(path) if isfile(join(path, f))]
faFiles = []
for i in faF:
    if "fasta.fai" not in i:
        faFiles.append(i)
df = pd.read_csv(file,sep=",")

housegenes = list(df.iloc[0:int(houseGeneEndRow),0])

df = pd.read_csv(file,sep=",")

# Finder samples
excludeCols = ['Sort column', 'Gene', 'contig_srt', 'contig_end', 'Target']
sampleList = []
discardedList = []

for col in df.columns:
  
    if col not in excludeCols:

        # Checking housegenes, first get total cov of them
        hgTot = 0
        for gene in housegenes:
            
            hg = df[df['Gene']==gene]
            hgTot+= int(hg[col])

        # Now finding outliers in housegene coverage
        notEnoughCov = 0
        for gene in housegenes:
            
            hg = df[df['Gene']==gene]
            geneCov = int(hg[col])

            # set min cov on housegene to not count as failed
            if geneCov < hgTot*0.1:
                notEnoughCov+=1
            
        # discard sample if too many low housegenes
        if notEnoughCov < 4:
            sampleList.append(col)
        else:
            discardedList.append(col)
discardedList

# Now remove housegenes from samples
ddf = df[~df['Gene'].isin(housegenes)]

ampl_df = pd.read_csv(amplPos,sep="\t", header = None)

# Getting amplicon id
amplList = []
for i in ddf['Gene']:
    amplList.append(i.split("_")[-1])
ddf['Target'] = amplList

ddf[['Type', 'Gene']] = ddf['Gene'].str.split('_', 1, expand=True)
# Removing Ampl from gene
# Getting amplicon id
GeneList = []
for i in ddf['Gene']:
    GeneList.append(i.split("AMPL")[~-1])
ddf['Gene'] = GeneList

ampl = pd.DataFrame()

ampl['Target'] = ampl_df[3]
ampl['contig_srt'] = ampl_df[1]
ampl['contig_end'] = ampl_df[2]

ddf = ddf.merge(ampl,on='Target')

# Look for deletions (possible integrations of the genome) by finding regions of very low coverage
hpvTypes = pd.unique(ddf['Type'])

for type in hpvTypes:
    filtered = ddf[ddf['Type'] == type]

delJson = {}

for sample in sampleList:
    for hpvtype in hpvTypes:
        filtered = ddf[ddf['Type'] == hpvtype]
        leng = len(filtered[sample])
        sum = filtered[sample].sum()
        
        # If sum of HPV type cover is less than n, do not call an integration
        if sum < 50:
            continue

        avg = sum / leng
        if avg > 50:
            cutoff = 50
        else:
            cutoff = avg*0.1

        #for sample in sampleList:#
        filtered = filtered.reset_index()
        lastIdx = 0
        for index, row in filtered.iterrows():
            
            if row[sample] < cutoff:
                # Check if last index contained low cov (potential deletion)
                if (lastIdx == index - 1) or (lastIdx == 0):

                    if 'startpos' not in locals():
                        startpos = row['contig_srt']
                    stoppos = row['contig_end']
                    lastIdx = index

                else:
                    # This is start of a new deletion
                    startpos = row['contig_srt']
                    lastrow = row
                    lastIdx = index


            else:
                # save found startpos and endpos
                if 'startpos' in locals():
                    if 'stoppos' not in locals():
                        stoppos = lastrow['contig_end']

                    # create key if it does not exist:
                    if sample in delJson.keys():
                        pass
                    else:    
                        delJson[sample] = {
                        'del' : []
                        }

                    delJson[sample]['del'].append([hpvtype,startpos,stoppos])

                    #save
                    del startpos, stoppos

integrationDF = pd.DataFrame(delJson)

filename = refSavePath+"/PossibleIntegrations.txt"

integrationDF.to_csv(filename)