import os, sys
import glob
import pandas as pd
from os import listdir
from os.path import isfile, join
from argparse import ArgumentParser

parser = ArgumentParser()

# define arguments. Argparse create help function from below help messages
parser.add_argument("-i", "--infile", required=True, help="Covarage matrix as csv file without qoutes, seperated by commas")
parser.add_argument("-r", "--ref", required=True, help="Location of all valid reference names that can be given")
parser.add_argument("-o","--outpath", required=True, help='Location of where to save output files')
parser.add_argument("-t","--typeTier", required=True, help='Tier of type to find. Maintype=1, subtype=2, etc.')

# to allow passing lists, include the nargs='+'. The list must be given as ' ' seperated values, eg. --inlist a b c 
# to set a defualt, use default='', to set another variable, use dest='myargs' e.g.
# This will then be callable as args.myarg. To set a type, use type=int e.g.

args = parser.parse_args()

# Housegenes to check and remove
housegenes = ['PABPN1','SRSF3','PPIE','RAB1B','BTF3']

#Covarage matrix as csv file without qoutes, seperated by commas
file = args.infile

# Location of all valid reference names that can be given
path = args.ref

# Location of where to save output files
refSavePath = args.outpath
# Remove all current inputRefs to avoid appending to them
if len(refSavePath) > 0:
    for f in glob.glob(os.path.join(refSavePath,'*.txt')):
        os.remove(f)


faF = [f for f in listdir(path) if isfile(join(path, f))]
faFiles = []
for i in faF:
    if "fasta.fai" not in i:
        faFiles.append(i)
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
# Split i samples og find hvilke HPV typer der passer
hpvTypes = pd.unique(ddf['Gene'])
# For sample, find HPV types in them
covJson = {}
for sample in sampleList:
    covJson.update({sample : {
    }})

    for hpv in hpvTypes:
        filtered = ddf[ddf['Gene'] == hpv]
        
        # Tilsætter sample, hpvtype og sum af dækning for regioner i hpvtype
        covJson[sample].update({
            hpv : filtered[sample].sum()
            })
# Tæl sum af dækning over alle regioner for sample, derefter angiv fraktioner

ratioJson = {}


for sample in covJson.keys():
    # Tæl total cov over gener for sample
    sampleCovTotal = 0
    for hpvtype, cov in covJson[sample].items():
        sampleCovTotal+=cov
    # Har nu total, og kan derfor loope igennem igen og lave ratios
    # Sætter sample name:
    ratioJson.update({sample : {
        }})
    
    for hpvtype, cov in covJson[sample].items():
        

        ratioJson[sample].update({
            hpvtype : cov/sampleCovTotal*100
            })
confirmedJson = {}


for sample in ratioJson.keys():
    
    confirmedJson.update({sample : {
        }})
    
    foundOneForSample = "no"
    highestRatio = 0
    for hpvtype, ratioCov in ratioJson[sample].items():
        
        # Define minimum ratio to confirm hpvtype in sample (1 to 100)
        if ratioCov > 20: 
            confirm = "yes" 
            foundOneForSample = "yes"
            if ratioCov > highestRatio:
                highestRatio = ratioCov
                highestRatioName = hpvtype
        else: 
            confirm = "no"

        confirmedJson[sample].update({
            hpvtype : confirm
            })

        # If no one meets the minimum ratioCov value, set the highest cov as best match
        if foundOneForSample == "no":
            for hpvtype, ratioCov in ratioJson[sample].items():
                if ratioCov > highestRatio:
                    highestRatio = ratioCov
                    highestRatioName = hpvtype
        
            for hpvtype, ratioCov in ratioJson[sample].items():
                if hpvtype == highestRatioName:
                    confirmedJson[sample].update({
                        hpvtype : "yes"
                        })



# For each sample output confirmed hpv types to file
multiFa = []
for sample in confirmedJson.keys():
    for hpvtype, val in confirmedJson[sample].items():
        if val == "yes":

            # convert name to fit fasta files:
            ht = hpvtype.split('_')[1]

            # search list of fasta files for id
            ref = [i for i in faFiles if ht in i]
            # Checking that there is only one match, then saving without .fasta in name
            if len(ref) > 1:
                multiFa.append([sample, ref])
            else:
                ref = ref[0].replace(".fasta","")
                # Append to file
                filename = refSavePath+"/"+sample+"T"+typeTier+".txt"
                f = open(filename, "a")
                f.write(ref+"\n")
                f.close()

# Printing fastas found with multiple hits. Fix these by removing one or more of the files
if len(multiFa) > 0:
    # Append to file
    filename = refSavePath+"/"+sample+"_foundMultiRefsFor.txt"
    f = open(filename, "a")
    f.write(ref+"\n")
    f.close()
