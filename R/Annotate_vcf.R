
# Predict amino acid changes from vcf, gff3 and fasta files

library(tidyverse)
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures", force = TRUE)
library(VariantAnnotation)
library(GenomicFeatures)

installed.packages()[, c("Package", "LibPath")]

##### EDIT ######
Fastqname <- "Pt_44.IonXpress_090_run"
Refname <- "K02718.1" # Uden extension
vcfname <- "K02718.1.align.vcf"
vcfout <- "K02718.1_anno.vcf" 
gffname <- "K02718.1.gff3"
#################

Ref <- paste("/home/pato/Skrivebord/HPV16_projekt/References/", Refname, ".fasta", sep = "")
faf <- open(FaFile(Ref))
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Sigma_run/", Fastqname, "/sigma_alignments_output/",Refname,"/", sep="")
vcffile <- paste(Folder, vcfname, sep ="")
gfffile <- paste(Folder, gffname, sep ="")
c_vcf <- readVcf(vcffile) # læser som collapsed vcf
gfffile

Gene_anno <- makeTxDbFromGFF(gfffile) # Laver TxDb objekt fra gff3 eller gtf fil. # , circ_seqs = "K02718.1"
# LIG MÆRKE TIL ADVARSELSBESKEDER

view(predictCoding(c_vcf, Gene_anno, seqSource = faf))

vcf_anno <- predictCoding(c_vcf, Gene_anno, seqSource = faf)
view(vcf_anno[,1])

vcf_annodf <- vcf_anno[,1]
view(vcf_annodf[1,1[1]])


#?predictCoding











  