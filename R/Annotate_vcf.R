
# Predict amino acid changes from vcf, gff3 and fasta files

#install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures", force = TRUE)
#BiocManager::install("Repitools")
library(tidyverse)
library(VariantAnnotation)
library(GenomicFeatures)
library(Repitools)
library(dplyr)

##### EDIT ######
Fastqname <- "Pt_33_RNA.IonXpress_087"
Refname <- "K02718.1" # Uden extension
Runname <- "Pt_33_RNA.IonXpress_087_run"  
vcfname <- paste(Refname,".vcf", sep ="")
gffname <- paste(Refname,"_E1fix.gff3", sep ="")
vcfout <- paste(Refname,"_anno.vcf", sep ="")
#################

Ref <- paste("/home/pato/Skrivebord/HPV16_projekt/References/", Refname, ".fasta", sep = "")
faf <- open(FaFile(Ref))
Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Results/", Runname,"/",Refname,"/", sep="")
GFFFolder <- paste("/home/pato/Skrivebord/HPV16_projekt/References/GFFfiles/", sep="")
vcffile <- paste(Folder, vcfname, sep ="")
gfffile <- paste(GFFFolder, gffname, sep ="")
c_vcf <- readVcf(vcffile) # læser som collapsed vcf

Gene_anno <- makeTxDbFromGFF(gfffile) # Laver TxDb objekt fra gff3 eller gtf fil. # , circ_seqs = "K02718.1"
# LIG MÆRKE TIL ADVARSELSBESKEDER

vcf_anno <- predictCoding(c_vcf, Gene_anno, seqSource = faf)
Nuc_start <- vcf_anno@ranges@start # Collecting nuc position
vcf_anno_df_2 <- vcf_anno@elementMetadata # Collecting elementdata
vcf_anno_df <- cbind(Nuc_start, vcf_anno_df_2)
view(vcf_anno_df)

anno_df <- as.data.frame(vcf_anno_df)
#anno_df[anno_df==""] <- NOTAVAIL

# Formatting protein changes
AAchange <- anno_df[c("REFAA","PROTEINLOC", "VARAA")] # "REFAA","PROTEIONLOC","VARAA"
AA <- cbind(p = "p.", AAchange)
rdyAA <- as.data.frame(paste(AA$p, AA$REFAA, AA$PROTEINLOC, AA$VARAA, sep = ""), col.names = aa)

# Formatting nucleotide changes
Nuchange <- vcf_anno_df[c("Nuc_start","REF", "ALT")] # "REFAA","PROTEIONLOC","VARAA"
Nuc <- cbind(Nuchange, a = ">")
Nuc <- cbind(c = "c.", Nuc)
rdyNuc <- as.data.frame(paste(Nuc$c, Nuc$Nuc_start, Nuc$REF, Nuc$a, Nuc$ALT@unlistData, sep = ""))
colnames(rdyNuc) <- "NuChange"

AnnoRdy <- cbind(Fastqname,rdyNuc,rdyAA,anno_df$GENEID)
colnames(AnnoRdy) <- c("Fastqname","Nuc","AA","GeneID")

view(AnnoRdy)

write.table(AnnoRdy, file = paste(Folder,"Annotation_", Fastqname,"_", Refname, ".txt", sep = ""), row.names = F,col.names = F, quote = F)






  
  