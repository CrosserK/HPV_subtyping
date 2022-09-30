
# Henter hver kolonne fra VCF filer som er relevant og ligger ud i SV_analysis/{SVCaller}

import math
import vcf
import pandas as pd
import time
import numpy as np

# tidsudløst, i sekunder:
#time.sleep(32400)


ReferenceList = ["deletions_d_12_100", "insertions_d_12_100", "duplications_d_12_100"]
MeanLengthList = ["15000","30000"]
TotalBasesList = ["500000000","1000000000","2000000000"]
AlignmentList = ["Minimap2","NGMLR"]
SVCallerList = ["SVIM"]
threads=12

                  # Finder filnavn
                  stdev = 5000 #int(meanLength) / 2
                  stdev = math.floor(stdev)  # Sørger for at tallet er uden decimaler

                  MainFolder = "/home/pato/Sekventerings_data/Kasper_speciale"
                  BamPath = MainFolder + "/Alignments/" + aligner
                  Fastqname = "chr17_" + refType + "_" + str(total_bases[:-6]) + "MbpData_" + str(
                      meanLength) + "Mean_" + str(
                      stdev) + "stdev.fastq.gz"
                  SortedBamFile = Fastqname[:-9] + ".Sorted." + aligner + ".bam"
                  ReferenceFolder = MainFolder + "/References"

                  # Finder korrekt vcf filnavn
                  VCF_outFolder = MainFolder + "/Alignments/" + aligner + "/SV_algorithm_outputs/" + SVCaller
                  VCF_name = "ReadForCall" + str(NumbOfReadsForCall) + "_MinLen" + str(
                      MinLengthForSV) + "_MinQual" + str(
                      Minmapping_qual) + "_" + SortedBamFile[:-4] + "_" + SVCaller + ".vcf"
                  VCF_out = VCF_outFolder + "/" + VCF_name
                  SortedVCF_out = MainFolder + "/Alignments/" + aligner + "/SV_algorithm_outputs/" + SVCaller + "/" + VCF_name[:-4] + "_sorted.vcf"

                  # Fixer navn hvis de ikke er sorted Sniffles eller NanoSV filer
                  if SVCaller == "Sniffles" or "NanoSV" or "cuteSV":
                      SortedVCF_out = SortedVCF_out[:-11] + ".vcf"

                  # Fixer navn hvis det er fixede SVIM filer
                  if SVCaller == "SVIM":
                      SortedVCF_out = SortedVCF_out[:-4] + "_sorted_fixed.vcf"



                  # Reading sv_subtype, svlen and startpos

                  try:
                    vcf_reader = vcf.Reader(open(SortedVCF_out, 'r'))
                  except:
                      print("Didn't find file" + SortedVCF_out)
                      fails = fails+1

                      # Printer fil med 0 0 0 0
                      """
                      data = {'first_set_of_numbers': [0],
                              'second_set_of_numbers': [0],
                              'third_set_of_numbers': [0],
                              'fourth_set_of_numbers': [0]
                              }
                      df = pd.DataFrame(data, columns=['first_set_of_numbers', 'second_set_of_numbers',
                                                       'third_set_of_numbers','fourth_set_of_numbers'])
                    
                      outfilename = SVCaller + "/" + VCF_name[:-4] + "_sorted.txt"
                      outfile = MainFolder + "/SV_analysis/" + outfilename
                      df.to_csv(outfile, index=False, sep="\t", header=False)
                      """
                  else:
                      svtypelist = []
                      svlenlist = []
                      svstartpos = []
                      svendpos = []
                      Stramd = []

                      # Checker om SVLEN er på listeform som fra NanoSV:
                      if SVCaller == "NanoSV":
                          for record in vcf_reader:
                              svtypelist.append(record.var_subtype)
                              svlenlist.append(abs(record.INFO["SVLEN"][0]))
                              svstartpos.append(record.POS)
                              svendpos.append(record.POS + abs(record.INFO["SVLEN"][0]))
                      else:
                          for record in vcf_reader:
                              svtypelist.append(record.var_subtype)
                              svlenlist.append(abs(record.INFO["SVLEN"]))
                              svstartpos.append(record.POS)
                              svendpos.append(record.POS + abs(record.INFO["SVLEN"]))

                      # Making dataframe:
                      SVdf = {'SV_type': svtypelist,
                              'SV_len': svlenlist,
                              'Startpos': svstartpos,
                              'Endpos' : svendpos
                              }

                      df = pd.DataFrame(SVdf,
                                        columns=['SV_type',
                                                 'SV_len',
                                                 'Startpos',
                                                 'Endpos'])  # I columns skal stå de samme navne som angivet i SVdf
                      df.sort_values(by='Startpos', inplace=True)
                      #print("Using:" + SortedVCF_out)
                      outfilename = SVCaller + "/" + VCF_name[:-4] + "_sorted.txt"
                      outfile = MainFolder + "/SV_analysis/" + outfilename
                      df.to_csv(outfile, index=False, sep="\t", header=False)
                      #print("Created: " + outfilename)

print("Could not be found: " + str(fails))
