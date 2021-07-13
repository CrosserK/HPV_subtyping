
# THIS SCRIP:
# Edit reference ID name to correct version after bp_genbank2gff processing

# Find all positions of char in string and correct VCF positions
library(Biostrings)
library(tidyverse)
### EDIT:
Fastqname <- "Pt_44.IonXpress_090_run"
gffname <- "K02718.1.gff3"
gffout <- "K02718.1_fixed.gff3" 
findtxt <- "K02718"
replacetxt <- "K02718.1"
#############
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Fastqname, "/sigma_alignments_output/K02718.1/", sep="") 
gff <- readLines(paste(Folder, gffname, sep = ""))
fixedgff  <- sub(pattern = findtxt, replace = replacetxt, x = gff)
writeLines(fixedgff, paste(Folder, gffout, sept = ""))
# Hvis der ikke sker noget når den trækkes ind i IGV, så omdøb den manuelt i mappen. 