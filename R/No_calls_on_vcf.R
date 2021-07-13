#!/usr/bin/env Rscript

# Køres efter Sitecoverage.sh for at lave no calls

library(tidyverse)
library(dplyr)

##TEST
# MainF <- "/home/pato/Skrivebord/HPV16_projekt"
# SaveDir <- "/home/pato/Skrivebord/HPV16_projekt/Annotation_results"
# SuperRunName <- "Karoline_run_covGenTest_0859_07072021"
# MultiFQfile <- paste("FASTQfiles_", SuperRunName, sep = "")
#####

MainF <- as.character(commandArgs(TRUE)[1])
SaveDir <- as.character(commandArgs(TRUE)[2])
SuperRunName <- as.character(commandArgs(TRUE)[3])
MultiFQfile <- as.character(commandArgs(TRUE)[4])

# Finder hver reference der er lavet data på skal der nu tjekkes om fastqfiler som er blevet alignet til reference også har dækning på position med snv
Refs <- list.files(path = SaveDir, pattern = paste("FASTQfiles_", SuperRunName, "_Nuc_change_coords_", sep =""))
Refs <- str_remove(Refs, paste("FASTQfiles_", SuperRunName, "_Nuc_change_coords_", sep =""))
Refs <- str_remove(Refs, ".bed")
# Har nu cuttet ned til referencer
# Loader ForNoCallScript fil
NoCalls <- read.table(paste(SaveDir,"/","ForNoCallScript_",MultiFQfile,"_Nuc_change_coords.txt", sep = ""), header = T)
# Laver nu no calls til hver reference
for(Refname in Refs){
  
  #TEST
  #Refname <- "K02718.1_revised"
  Fredf <- NoCalls[NoCalls$Reference == Refname,]
  
  MultiFastqListFile <- paste(MainF,"/", MultiFQfile, ".txt", sep = "")
  MultiFastqList <- read.table(MultiFastqListFile)
  MultiFastqList <- as.list(MultiFastqList)
  MultiFastqList <- unlist(MultiFastqList)
  
  lengthofList <- length(MultiFastqList)
  
  ##### Henter og tilsætter NO CALLS ######
  #FreqDF <- separate(Fredf, col = 2, sep = "-", into = c("Nuc", "-"))
  FreqDF <- Fredf
  counter <- 0 
  nocalldf <- data.frame()
  FreqDF <- cbind(FreqDF, No_Calls = 0)
  FreqDF$No_Calls <- as.integer(FreqDF$No_Calls)
  
  for(Fastqname in MultiFastqList){
    
    #TEST
    #Fastqname <- "pt_130.IonXpress_006"
    #####
   
    # Tjekker om fastqfil er blevet alignet til nuværende reference:
    if(file.exists(paste(MainF,"/Results/",SuperRunName,"/", Fastqname,"/",Refname,"/ResultFiles/SNP_cov.txt", sep="")) == FALSE){
      next
    } else {
      print("yes")
    }
    
    Runname <- paste(Fastqname,"_run",sep="")
    Folder <- paste(MainF,"/Results/",SuperRunName,"/", Fastqname,"/",Refname,"/ResultFiles/", sep="")

    
    # Hopper til næste iteration i loop, hvis ingen SNP_cov fil (Da der så ikke er blevet alignet til nuværende reference med nuværende fastqfil)
    SNP_cov_File <- try(read.table(paste(Folder,"SNP_cov.txt", sep ="")))
    if("try-error" %in% class(SNP_cov_File)){
      next
    }
    
    for(j in 1:nrow(FreqDF)){
      for(n in 1:nrow(SNP_cov_File)){
        if((FreqDF[j,"NucPos"] == SNP_cov_File[n,2]) & (SNP_cov_File[n,3] < 5) ){
          FreqDF[j,"No_Calls"] <- FreqDF[j,"No_Calls"] + 1
        }
      }
    }
    
    counter <- counter + 1
    
    # Convenience progress tracker
    print(paste(Fastqname,"is done...",counter,"of",lengthofList, sep = " "))
    if(counter == lengthofList){
      print("Job's done!")
      
    }
  }
  
  # Printer annotation for nuværende reference der nu er testet på alle fastq filer
  
  Anno_freq <- FreqDF 
  
  # Samler tabel til pæn form  
  Anno_freq[,"NucChange"] <- as.data.frame(paste("c.",Anno_freq$NucPos, Anno_freq$NucChange, sep = ""))
  Anno_freq <- Anno_freq[,!names(Anno_freq) %in% "NucPos"]
  # Laver annotation info filer:
  write.table(Anno_freq, file = paste(SaveDir,"/","AnnotationFrequency_",MultiFQfile, Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F, sep = "\t")
  #file.remove(paste(SaveDir,"/","ForNoCallScript_",MultiFQfile,"_Nuc_change_coords_", Refname, ".txt", sep = ""))
}

fn <- paste(SaveDir, "/","ForNoCallScript_",MultiFQfile,"_Nuc_change_coords.txt", sep = "")

# fjerner NoCallFile fra annoscript
if (file.exists(fn)) {
  #Delete file if it exists
  file.remove(fn)
}



