
# Predict amino acid changes from vcf, gff3 and fasta files

#install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures", force = TRUE)
#BiocManager::install("Repitools")
library(tidyverse)
library(VariantAnnotation)
library(GenomicFeatures)
library(Repitools)

##### EDIT ######
Fastqname <- "Pt_43.IonXpress_089"
Refname <- "AF402678.1" # Uden extension
Runname <- "run_001"  
vcfname <- paste(Refname,".align.vcf", sep ="")
gffname <- paste(Refname,".gff3", sep ="")
vcfout <- paste(Refname,"_anno.vcf", sep ="") 
#################

Ref <- paste("/home/pato/Skrivebord/HPV16_projekt/References/", Refname, ".fasta", sep = "")
faf <- open(FaFile(Ref))
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Runname, "/sigma_alignments_output/",Refname,"/", sep="")
vcffile <- paste(Folder, vcfname, sep ="")
gfffile <- paste(Folder, gffname, sep ="")
c_vcf <- readVcf(vcffile) # læser som collapsed vcf
gfffile

Gene_anno <- makeTxDbFromGFF(gfffile) # Laver TxDb objekt fra gff3 eller gtf fil. # , circ_seqs = "K02718.1"
# LIG MÆRKE TIL ADVARSELSBESKEDER

view(predictCoding(c_vcf, Gene_anno, seqSource = faf))

vcf_anno <- predictCoding(c_vcf, Gene_anno, seqSource = faf)



# Formaterer
vcf_anno_df <- annoGR2DF(vcf_anno)
AAchange <- select(vcf_anno_df, REFAA, PROTEIONLOC, VARAA) # "REFAA","PROTEIONLOC","VARAA"
NucChange <- select(vcf_anno_df, CDSLOC.start, REF, ALT)








  