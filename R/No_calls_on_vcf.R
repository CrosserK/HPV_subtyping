#!/usr/bin/env Rscript

# Køres efter Sitecoverage.sh for at lave no calls

library(tidyverse)
library(dplyr)

SaveDir <- as.character(commandArgs(TRUE)[1])
Refname <- as.character(commandArgs(TRUE)[2])
SuperRunName <- as.character(commandArgs(TRUE)[3])
MultiFQfile <- as.character(commandArgs(TRUE)[4])

#SaveDir <- "/home/pato/Skrivebord/HPV16_projekt/Annotation_results/"
#Refname <- "K02718.1_PaVE"
#GFFname <- "K02718.1_revised"
#SuperRunName <- "Exome_50_120_ampliconcalls_PaVE_revised"
#MultiFQfile <- "FASTQfiles_Exome_50_120_ampliconcalls_PaVE"

Fredf <- read.table(paste(SaveDir,"ForNoCallScript_",MultiFQfile,"_Nuc_change_coords_", Refname, ".txt", sep = ""), header = T)

MultiFastqListFile <- paste("/home/pato/Skrivebord/HPV16_projekt/", MultiFQfile, ".txt", sep = "")
MultiFastqList <- read.table(MultiFastqListFile)
MultiFastqList <- as.list(MultiFastqList)
MultiFastqList <- unlist(MultiFastqList) 

lengthofList <- length(MultiFastqList)

##### Henter og tilsætter NO CALLS ######
FreqDF <- separate(Fredf, col = 2, sep = "-", into = c("Nuc", "-"))
counter <- 0 
nocalldf <- data.frame()
FreqDF <- cbind(FreqDF, No_Calls = 0)
FreqDF$No_Calls <- as.integer(FreqDF$No_Calls)
for(i in MultiFastqList){
  
  Fastqname <- i
  Runname <- paste(i,"_run",sep="")
  Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Results/",SuperRunName,"/", Runname,"/",Refname,"/ResultFiles/", sep="")
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
write.table(Anno_freq, file = paste(SaveDir,MultiFQfile,"_AnnotationFrequency_", Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F)

