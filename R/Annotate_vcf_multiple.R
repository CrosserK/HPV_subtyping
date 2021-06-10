
# Predict amino acid changes from vcf, gff3 and fasta files, this version will output format for multiple fastqfiles
# in a convinient way

#install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures", force = TRUE)
#BiocManager::install("Repitools")
library(tidyverse)
library(VariantAnnotation)
library(GenomicFeatures)
library(Repitools)
library(dplyr)

########## EDIT ##############
SaveDir <- "/home/pato/Skrivebord/HPV16_projekt/Annotation_results/"
Refname <- "K02718.1"
SuperRunName <- "Exome_50_120"
MultiFQfile <- "FASTQfiles_Exome_50_120"
##############################

MultiFastqListFile <- paste("/home/pato/Skrivebord/HPV16_projekt/", MultiFQfile, ".txt", sep = "")
MultiFastqList <- read.table(MultiFastqListFile)
MultiFastqList <- as.list(MultiFastqList)
MultiFastqList <- unlist(MultiFastqList) # Unlister for at for loop kan læse korrekt


vcfname <- paste(Refname,"_filtered.filtEx_headerfix.vcf", sep ="") # _filtered.filtEx_headerfix
gffname <- paste(Refname,"_E1fix.gff3", sep ="")
vcfout <- paste(Refname,"_anno.vcf", sep ="") 


cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
  rbind(x, matrix(,n-nrow(x), ncol(x))))) 
}
GFFFolder <- paste("/home/pato/Skrivebord/HPV16_projekt/References/GFFfiles/", sep="")
gfffile <- paste(GFFFolder, gffname, sep ="")
Gene_anno <- makeTxDbFromGFF(gfffile) # Laver TxDb objekt fra gff3 eller gtf fil. # , circ_seqs = "K02718.1"
# LIG MÆRKE TIL ADVARSELSBESKEDER
newdf <- data.frame()
newdf_nuc <- data.frame()
newdf_nuc_coord <- data.frame()
lengthofList <- length(MultiFastqList)
counter <- 0
noelement <- 0
novcfcounter <- 0
for(i in MultiFastqList){
  ###TEST
  #i <- "Pt_79_RNA.IonXpress_073"
  #####
  
  Fastqname <- i
  Runname <- paste(i,"_run",sep="")
  
  Ref <- paste("/home/pato/Skrivebord/HPV16_projekt/References/", Refname, ".fasta", sep = "")
  faf <- open(FaFile(Ref))
  Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Results/",SuperRunName,"/", Runname,"/",Refname,"/", sep="")
  counter <- counter + 1
  # Error check for om der er en filtreret vcf fil tilgængelig
  vcffile <- paste(Folder, vcfname, sep ="")
  
  c_vcf <- try(readVcf(vcffile)) # læser som collapsed vcf
  if("try-error" %in% class(c_vcf)){
    novcfcounter <- novcfcounter + 1
    next
  }

  # Error check som fortsætter loopet hvis der ingen varianter i vcf fil er
  options(warn=2)
  no_elem <- try(predictCoding(c_vcf, Gene_anno, seqSource = faf))
  if("try-error" %in% class(no_elem)){
    noelement <- noelement + 1
    next
  }
  options(warn=1)  

  vcf_anno <- predictCoding(c_vcf, Gene_anno, seqSource = faf)
  Nuc_start <- vcf_anno@ranges@start # Collecting nuc position
  vcf_anno_df_2 <- vcf_anno@elementMetadata # Collecting elementdata
  vcf_anno_df <- cbind(Nuc_start, vcf_anno_df_2)

  anno_df <- as.data.frame(vcf_anno_df)
  #anno_df[anno_df==""]
  
  # Formatting protein changes
  AAchange <- anno_df[c("REFAA","PROTEINLOC", "VARAA","GENEID")] # "REFAA","PROTEIONLOC","VARAA"
  AA <- cbind(p = "p.", AAchange)
  rdyAA <- as.data.frame(paste(AA$p, AA$REFAA, AA$PROTEINLOC, AA$VARAA,"_", AA$GENEID, sep = ""), col.names = aa)

  # Formatting nucleotide changes
  Nuchange <- vcf_anno_df[c("Nuc_start","REF", "ALT")] # "REFAA","PROTEIONLOC","VARAA"
  Nuc <- cbind(Nuchange, a = ">")
  #Nuc <- cbind(c = "c.", Nuc)
  Nuc <- cbind(sepera = ",", Nuc)
  rdyNuc <- as.data.frame(paste(Nuc$Nuc_start,Nuc$sepera, Nuc$REF, Nuc$a, Nuc$ALT@unlistData, sep = ""))
  colnames(rdyNuc) <- "Mutation"
  rm(AnnoRdy) 
  #Merger AA og nuc info til 1 kolonne
  AnnoRdy <- data.frame(matrix(ncol = 1, nrow = nrow(rdyAA)))
  for(j in 1:nrow(rdyAA)){
    AnnoRdy[j,] <- paste(rdyNuc[j,],rdyAA[j,],sep = "-") 
  }
    
  colnames(AnnoRdy) <- Fastqname
  
  newdf <- cbind.fill(newdf, AnnoRdy)
  
  #newdf_nuc <- cbind.fill(newdf_nuc, rdyNuc)
  
  #newdf_nuc_coord <- cbind.fill(newdf_nuc_coord,Nuc$Nuc_start)
  
  # Conveniece for tracking progress
  print(paste(i,"is done...",counter,"of",lengthofList, sep = " "))
  
}

print("Job's done!")
print(paste("Number of noelements found: ", noelement))
print(paste("Number of nofiltered vcf: ", novcfcounter))

# Find each unique instance of change:
Frequency <- table(newdf)
Fredf <- as.data.frame(Frequency)
# Splitting data table into separate columns
Fredf <- separate(Fredf,col = 1, sep = "_", into = c("Change", "GENEID"))
Fredf <- separate(Fredf, col = 1, sep = ",", into = c("NucPos","Nuc"))
Fredf$NucPos <- as.integer(Fredf$NucPos)
Fredf <- Fredf[order(Fredf$NucPos),]
# Laver en bed fil opsætning med nuc change coords
NucChangePos <- Fredf$NucPos
df <- as.vector(as.matrix(NucChangePos))
# df <- unique(df) # Finder alle unikke værdier
df <- data.frame(Reference = Refname, Pos = NucChangePos-1, PosEnd = NucChangePos)
write.table(df, file = paste(SaveDir,MultiFQfile,"_Nuc_change_coords_", Refname, ".bed", sep = ""), row.names = F,col.names = F, quote = F, sep = "\t")





# Kør nu Sitecoverage.sh





##### Henter og tilsætter NO CALLS ######
FreqDF <- separate(Fredf, col = 2, sep = "-", into = c("Nuc", "-"))
counter <- 0 
nocalldf <- data.frame()
FreqDF <- cbind(FreqDF, No_Calls = 0)
FreqDF$No_Calls <- as.integer(FreqDF$No_Calls)
for(i in MultiFastqList){
  
  Fastqname <- i
  Runname <- paste(i,"_run",sep="")
  Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Results/",SuperRunName,"/", Runname,"/",Refname,"/", sep="")
  SNP_cov_File <- read.table(paste(Folder,"SNP_cov.txt", sep =""))
  
  for(j in 1:nrow(FreqDF)){
    for(n in 1:nrow(SNP_cov_File)){
      if((FreqDF[j,1] == SNP_cov_File[n,2]) & (SNP_cov_File[n,3] < 5) ){
        FreqDF[j,6] <- FreqDF[j,6] + 1
      }
    }
  }
  
  counter <- counter + 1

  # Convenience progress tracker
  print(paste(i,"is done...",counter,"of",lengthofList, sep = " "))
  if(counter == lengthofList){
    print("Job's done!")

  }
}

# Printer annotation
  
Anno_freq <- FreqDF 
Anno_freq[,1] <- as.data.frame(paste("c.",Anno_freq$NucPos, Anno_freq$Nuc, sep = ""))
Anno_freq <- Anno_freq[,-2]
colnames(Anno_freq) <- c("Base", "AA", "GENEID", "Freq", "No_call")
# Laver annotation info filer:
write.table(newdf, file = paste(SaveDir,MultiFQfile,"_Annotation_", Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F, sep = "\t")
write.table(Anno_freq, file = paste(SaveDir,MultiFQfile,"_AnnotationFrequency_", Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F)



