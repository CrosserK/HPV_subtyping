# Graph vcf stats

library(tidyverse)

################
Fastqname <- "Pt_79_RNA.IonXpress_092_run"
Refname <- "AF402678.1" # Uden extension
Runname <- "Pt_79_RNA.IonXpress_092_run"
# Choose stats to graph, currently only plots unfiltered files if run!
filename <- "depth.txt"
filteredfilename <- "depth_filt.txt"
#################

par(mfrow=c(2,2)) 

# Getting vcfstats
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Runname, "/sigma_alignments_output/",Refname,"/", sep="")

# Getting files
vcf <- paste(Folder,"/VCFstats","/", filename, sep ="")
vcf_filt <- paste(Folder,"/VCFstats","/", filteredfilename, sep ="")
vcf <- read.table(vcf)
vcf_filt <- read.table(vcf_filt)

# Plotting
plot(density(vcf$V1), main = "Depth")

filename <- "FS.txt"
filteredfilename <- "FS_filt.txt"
# Getting vcfstats
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Runname, "/sigma_alignments_output/",Refname,"/", sep="")

# Getting files
vcf <- paste(Folder,"/VCFstats","/", filename, sep ="")
vcf_filt <- paste(Folder,"/VCFstats","/", filteredfilename, sep ="")
vcf <- read.table(vcf)
vcf_filt <- read.table(vcf_filt)

# Plotting
plot(density(vcf$V1), main = "FS")

filename <- "MQ.txt"
filteredfilename <- "MQ_filt.txt"
# Getting vcfstats
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Runname, "/sigma_alignments_output/",Refname,"/", sep="")

# Getting files
vcf <- paste(Folder,"/VCFstats","/", filename, sep ="")
vcf_filt <- paste(Folder,"/VCFstats","/", filteredfilename, sep ="")
vcf <- read.table(vcf)
vcf_filt <- read.table(vcf_filt)

# Plotting
plot(density(vcf$V1), main = "MQ")


filename <- "QD.txt"
filteredfilename <- "QD_filt.txt"
# Getting vcfstats
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Runname, "/sigma_alignments_output/",Refname,"/", sep="")

# Getting files
vcf <- paste(Folder,"/VCFstats","/", filename, sep ="")
vcf_filt <- paste(Folder,"/VCFstats","/", filteredfilename, sep ="")
vcf <- read.table(vcf)
vcf_filt <- read.table(vcf_filt)

# Plotting
plot(density(vcf$V1), main = "QD")
