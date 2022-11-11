import os, sys
import glob
import pandas as pd
from os import listdir
from os.path import isfile, join
from argparse import ArgumentParser

parser = ArgumentParser()

# define arguments. Argparse create help function from below help messages
parser.add_argument("-i", "--infile", required=True, help="Coverage matrix as csv file without qoutes, seperated by commas")
parser.add_argument("-r", "--ref", required=True, help="Location of all valid references that can be linked to")
parser.add_argument("-f", "--fastqFiles", required=True, help="Location of all fastqfiles")
parser.add_argument("-o","--outpath", required=True, help='Location of where to save output files')
parser.add_argument("-t","--typeTier", required=True, help='Tier of type to find. Maintype=1, subtype=2, etc.')
parser.add_argument("-g","--houseGeneEndRow", required=True, help='Last row containing housegenes. They should be from row 1 to this row')
parser.add_argument("-l","--logFile", required=True, help='Location of log file')

# to allow passing lists, include the nargs='+'. The list must be given as ' ' seperated values, eg. --inlist a b c 
# to set a defualt, use default='', to set another variable, use dest='myargs' e.g.
# This will then be callable as args.myarg. To set a type, use type=int e.g.

args = parser.parse_args()

# Housegenes to check and remove
#housegenes = ['PABPN1','SRSF3','PPIE','RAB1B','BTF3']
#housegenes = ['chr1_PPIE_AMPL4425884','chr5_BTF3_AMPL4425932','chr6_SRSF3_AMPL4425918','chr11_RAB1B_AMPL4426066','chr14_PABPN1_AMPL4425941']


#Covarage matrix as csv file without quotes, seperated by commas
file = args.infile
typeTier = args.typeTier
fastqFiles = args.fastqFiles

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

housegenes = list(df.iloc[0:int(args.houseGeneEndRow),0])
print("Housegenes: ",housegenes)
f = open(args.logFile, "a")
f.write("Housegenes set in covMatrix: "+str(housegenes)+"\n")
f.close()

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

# Split type-gene column 
ddf[['Type', 'Gene']] = ddf['Gene'].str.split('_', 1, expand=True)

# Find unique hpv types in matrix
hpvTypes = pd.unique(ddf['Type'])

# For sample, find HPV types in them
covJson = {}
for sample in sampleList:
    covJson.update({sample : {
    }})

    for hpv in hpvTypes:
        filtered = ddf[ddf['Type'] == hpv]
        
        # Checking how many amplicons have more than n reads
        ampliconsWithReads = 0
        for i in filtered[sample]:
            # define minimum reads to count amplicon as covered
            if i > 1: 
                ampliconsWithReads+=1
        # Define minimum number of amplicons to cover
        if ampliconsWithReads < 2:
            covJson[sample].update({
                hpv : 0
                })
        else:
            # Tilsætter sample, hpvtype og sum af dækning for regioner i hpvtype
            covJson[sample].update({
                hpv : filtered[sample].sum()
                })

# Tæl sum af dækning over alle regioner for sample, derefter angiv fraktioner
ratioJson = {}
for sample in covJson.keys():
    # Tæl total cov over gener for sample
    sampleCovTotal = 0
    ampliconsWithReads = 0
    for hpvtype, cov in covJson[sample].items():
        sampleCovTotal+=cov
    # Har nu total, og kan derfor loope igennem igen og lave ratios
    # Sætter sample name:
    ratioJson.update({sample : {
        }})
    
    for hpvtype, cov in covJson[sample].items():

        # Avoid division by 0:
        if cov == 0:
            ratioJson[sample].update({
            hpvtype : 0
            })
        else:
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

# For each sample output confirmed hpv types to file
multiFa = []
for sample in confirmedJson.keys():

    # Because the fastq are usually with a prefix compared to in the covMatrix, find the fastq name from the matrix name
    for file in glob.glob(fastqFiles+"/*"+sample+"*"):
        fastqname = file
        fastqname = os.path.basename(fastqname)
        fastqname = fastqname.replace(".fastq","")
    if "fastqname" not in locals():
        fastqname = sample + "_FASTQ_NOT_FOUND"

    file1 = open(refSavePath+"/TypeCallSummary_T1.txt", "a")
    L = [fastqname] 
    file1.writelines(L)
    file1.close()

    filename = refSavePath+"/"+fastqname+"/TypeCalls/"+fastqname+"_T"+str(typeTier)+"_SplitTo.txt"

    os.makedirs(refSavePath+"/"+fastqname+"/TypeCalls/", exist_ok = True)
    for hpvtype, val in confirmedJson[sample].items():
        if val == "yes":

            # Append _ to hpv name so that HPV161 does not match HPV16 in a search:
            if "_" in hpvtype:
                pass
            else:
                ht = hpvtype + "_" 

            # search list of fasta files for id
            ref = [i for i in faFiles if ht in i]
            # Checking that there is only one match, then saving without .fasta in name
            if len(ref) > 1:
                multiFa.append([sample, ref])
            else:
                # Checking that there is atleast one match
                if ref == []:
                    ref = "No_type_found"
                else:
                    ref = ref[0].replace(".fasta","")
                # Append to file
                f = open(filename, "a")
                f.write(ref+"\n")
                f.close()

                # Appending to summary
                file1 = open(refSavePath+"/TypeCallSummary_T1.txt", "a")
                L = ["\t", ref] 
                file1.writelines(L)
                file1.close()
                
    file1 = open(refSavePath+"/TypeCallSummary_T1.txt", "a")
    file1.writelines("\n")
    file1.close()
    del fastqname


# Printing fastas found with multiple hits. Fix these by removing one or more of the files
if len(multiFa) > 0:
    # Append to file
    filename = refSavePath+"/"+sample+"_foundMultiRefsFor.txt"
    f = open(filename, "a")
    f.write(ref+"\n")
    f.close()
