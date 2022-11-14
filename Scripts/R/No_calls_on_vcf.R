#!/usr/bin/env Rscript

# Køres efter Sitecoverage.sh for at lave no calls
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))

##TEST
# MainF <- "/home/pato/Skrivebord/HPV_subtyping"
# TopRunName <- "Reanalyze_user_S5-0184-347-HPV_genotypering_Sara_Opsaet_10_A_1459_11112022"
# MultiFQfile <- paste("/home/pato/Skrivebord/HPV_subtyping/FASTQ/FASTQfiles_",TopRunName,".txt",sep="")
# SaveDir <- paste("/home/pato/Skrivebord/HPV_subtyping/Results/",TopRunName,sep="")
# typeTier <- "2"
#####

MainF <- as.character(commandArgs(TRUE)[1])
SaveDir <- as.character(commandArgs(TRUE)[2])
TopRunName <- as.character(commandArgs(TRUE)[3])
MultiFQfile <- as.character(commandArgs(TRUE)[4])

# Finder hver reference der er lavet data på skal der nu tjekkes om fastqfiler som er blevet alignet til reference også har dækning på position med snv
Refs <- list.files(path = SaveDir, pattern = paste("VariantPositions_",TopRunName, ".bed", sep = ""))
snp_cov <- read.table(paste(SaveDir,"/",Refs,sep=""),sep="\t")
Refs <- unique(snp_cov[1])
#Reset index
rownames(Refs) <- NULL  
# Har nu cuttet ned til referencer
# Loader ForNoCallScript fil
NoCalls <- read.table(paste(SaveDir,"/","ForSummaryScript_",TopRunName,"_Nuc_change_coords.txt", sep = ""), header = T)
# Laver nu no calls til hver reference. +1 because the last one is with the NAs

for(i in 1:as.integer(length(Refs[[1]])+1)){
  # This is for getting NAs if on last iteration that is +1 of length
  if(i!=length(Refs[[1]])+1){
    Refname <- Refs[[1]][i]
    FreqDF <- NoCalls[NoCalls$Reference == Refname,]
    FreqDF <- na.omit(FreqDF)#FreqDF[FreqDF$GENEID != is.na(FreqDF$GENEID)]
  } else {
    Refname <- "isNA"
    FreqDF <- NoCalls[is.na(NoCalls$AAChange),]
  }
  
  
  #TEST
  #Refname <- "K02718.1_revised"
  #Refname <- "FN907962_1"
  

  # if(length(FreqDF)==0){
  #   next
  # }
  MultiFastqListFile <- MultiFQfile
  MultiFastqList <- read.table(MultiFastqListFile)
  MultiFastqList <- as.list(MultiFastqList)
  MultiFastqList <- unlist(MultiFastqList)
  
  lengthofList <- length(MultiFastqList)
  
  ##### Henter og tilsætter NO CALLS ######
  #FreqDF <- separate(FreqDF, col = 2, sep = "-", into = c("Nuc", "-"))
  counter <- 0 
  nocalldf <- data.frame()
  #FreqDF <- cbind(FreqDF, No_Calls = 0)
  #FreqDF$No_Calls <- as.integer(FreqDF$No_Calls)
  
  for(Fastqname in MultiFastqList){
    
    #TEST
    #Fastqname <- "pt_130.IonXpress_006"
  
    #####
    Fastqname <- tools::file_path_sans_ext(basename(Fastqname)) # Getting fastqname without path and ext
   
    Folder <- paste(MainF,"/Results/",TopRunName,"/", Fastqname,"/",Refname,"/ResultFiles/", sep="")
    
    # Hopper til næste iteration i loop, hvis ingen SNP_cov fil (Da der så ikke er blevet alignet til nuværende reference med nuværende fastqfil)
    #SNP_cov_File <- snp_cov
    #if("try-error" %in% class(SNP_cov_File)){
    #  next
    #}
    
    # for(j in 1:nrow(FreqDF)){
    #   for(n in 1:nrow(SNP_cov_File)){
    #     if((FreqDF[j,"NucPos"] == SNP_cov_File[n,2]) & (SNP_cov_File[n,3] < 5) ){
    #       FreqDF[j,"No_Calls"] <- FreqDF[j,"No_Calls"] + 1
    #     }
    #   }
    # }
    
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
  
  # Making table more convenient
  Anno_freq$GENEID[is.na(Anno_freq$GENEID)] <- "No suitable GFF file found"
  #Anno_freq$GENEID <- Anno_freq[is.na(Anno_freq$GENEID)] 
  Anno_freq$GENEID <- str_replace(Anno_freq$GENEID,"1,","")
  # Laver annotation info filer:
  write.table(Anno_freq, file = paste(SaveDir,"/","AnnotationFrequency_",TopRunName, Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F, sep = "\t")
  #file.remove(paste(SaveDir,"/","ForSummaryScript_",MultiFQfile,"_Nuc_change_coords_", Refname, ".txt", sep = ""))
}



# Cleaning up, by deleting no call file
fn <- paste(SaveDir, "/","ForSummaryScript_",TopRunName,"_Nuc_change_coords.txt", sep = "")
if (file.exists(fn)) {
 file.remove(fn)
}
print("Summaries done")

