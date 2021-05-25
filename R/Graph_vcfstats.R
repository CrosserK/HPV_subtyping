# Graph vcf stats

library(tidyverse)

################
Fastqname <- "Pt_44.IonXpress_090_run"
Refname <- "AF402678.1" # Uden extension
Runname <- "run_001"
# Choose stats to graph
filename <- "MQ.txt"
filteredfilename <- "MQ_filt.txt"
#################

# Getting vcfstats
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Runname, "/sigma_alignments_output/",Refname,"/", sep="")

# Getting files
vcf <- paste(Folder,"/VCFstats","/", filename, sep ="")
vcf_filt <- paste(Folder,"/VCFstats","/", filteredfilename, sep ="")
vcf <- read.table(vcf)
vcf_filt <- read.table(vcf_filt)

# Plotting
plot(vcf$V1)
plot(density(vcf$V1))



