import os, sys
import glob
import pandas as pd
from os import listdir
from os.path import isfile, join
from argparse import ArgumentParser

parser = ArgumentParser()

# define arguments. Argparse create help function from below help messages
parser.add_argument("-r", "--ref", required=True, help="Location of all valid references that can be linked to")
parser.add_argument("-o","--outpath", required=True, help='Location of where to save output files')
parser.add_argument("-g","--houseGeneEndRow", required=True, help='Last row containing housegenes. They should be from row 1 to this row')
parser.add_argument("-a","--amplPos", required=True, help='bed file with amplicon positions. See manual for format')

# to allow passing lists, include the nargs='+'. The list must be given as ' ' seperated values, eg. --inlist a b c 
# to set a defualt, use default='', to set another variable, use dest='myargs' e.g.
# This will then be callable as args.myarg. To set a type, use type=int e.g.

args = parser.parse_args()

# Location of all valid reference names that can be given
path = args.ref

# Location of where to save output files
outpath = args.outpath

# Find possible references to match
faF = [f for f in listdir(path) if isfile(join(path, f))]
faFiles = []
for i in faF:
    if "fasta.fai" not in i:
        faFiles.append(i)

ampl_df = pd.read_csv(args.amplPos,sep="\t", header = None)
ampl_df = ampl_df.iloc[int(args.houseGeneEndRow):,]
ampl_df[['Pool', 'Submitted_reg']] = ampl_df[5].str.split(';', 1, expand=True)

bedpool = []
for pool in ampl_df.Pool.unique():
    poolN = pool.split("=")[1]
    savename = "pool"+str(poolN)+".bed"
    #Filtering for pool
    filtered = ampl_df[ampl_df['Pool']==pool]
    # Cutting down to only hpvtype, start and end positions
    filtered = filtered.iloc[:,0:3]

    # Now create bed file with corrected hpvtype names for each pool
    corrected_hpvtypes = []
    for index, row in filtered.iterrows():
        hpvtype = row[0]
        contig_srt = row[1]
        contig_end = row[2]
        # Create dataframe with correct ID's
        # Append _ to hpv name so that HPV161 does not match HPV16 in a search:
        if "_" in hpvtype:
            ht = hpvtype.split("_")[1]
            pass
        else:
            ht = hpvtype + "_" 

        # search list of fasta files for id
        ref = [i for i in faFiles if ht in i]
        # Checking that there is only one match, then saving without .fasta in name
        if len(ref) > 1:
            corrected_hpvtypes.append(["MultipleRefsFound", contig_srt, contig_end])
        else:
            # Checking that there is atleast one match
            if ref == []:
                ref = "No_type_found"
            else:
                ref = ref[0].replace(".fasta","")
                corrected_hpvtypes.append([ref, contig_srt, contig_end])

    df = pd.DataFrame(corrected_hpvtypes)
    df.to_csv(outpath+"/"+savename, sep="\t", header = None, index = None)