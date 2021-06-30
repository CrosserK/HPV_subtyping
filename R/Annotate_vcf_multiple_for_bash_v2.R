#!/usr/bin/env Rscript

# Predict amino acid changes from vcf, gff3 and fasta files, this version will output format for multiple fastqfiles
# in a convenient way

#install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures", force = TRUE)
#BiocManager::install("Repitools")
library(GenomicFeatures)
library(tidyverse)
library(VariantAnnotation)
library(Repitools)
library(dplyr)

####TEST 
# MainF <- "/home/pato/Skrivebord/HPV16_projekt"
# SaveDir <- "/home/pato/Skrivebord/HPV16_projekt/Annotation_results"
# SuperRunName <- "Karoline_run_1007_30062021"
# MultiFQfile <- paste("FASTQfiles_", SuperRunName, sep = "")
######

# Henter variabler fra commandline argumenter
MainF <- as.character(commandArgs(TRUE)[1])
SaveDir <- as.character(commandArgs(TRUE)[2])
SuperRunName <- as.character(commandArgs(TRUE)[3])
MultiFQfile <- as.character(commandArgs(TRUE)[4])


MultiFastqListFile <- paste(MainF,"/", MultiFQfile, ".txt", sep = "")
MultiFastqList <- read.table(MultiFastqListFile)
MultiFastqList <- as.list(MultiFastqList)
MultiFastqList <- unlist(MultiFastqList) # Unlister for at for loop kan læse korrekt

# Laver funktion som gør at kolonner af forskellig længde kan sættes sammen, hvor kortere kolonner bliver fyldt op med "NA"
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
  rbind(x, matrix(,n-nrow(x), ncol(x))))) 
}
GFFFolder <- paste("/home/pato/Skrivebord/HPV16_projekt/References/GFFfiles/", sep="")

# LIG MÆRKE TIL ADVARSELSBESKEDER
newdf <- data.frame()
newdf_nuc <- data.frame()
newdf_nuc_coord <- data.frame()
lengthofList <- length(MultiFastqList)
counter <- 0
noelement <- 0
novcfcounter <- 0
for(Fastqname in MultiFastqList){
  ###TEST
  #Fastqname <- "Pt_79_RNA.IonXpress_092"
  #Fastqname <- "pt_45.IonXpress_067"
  #####
  
  Runname <- paste(Fastqname,"_run",sep="")
  # Finder hvilke referencer fastq filer er blevet alignet til 
  RevReferences <- read.table(paste(MainF,"/", "VirStrain_run/", SuperRunName,"/", Runname,"/","Revised_SubTypeCalls.txt", sep = ""))
  # Hopper til næste iteration (næste fastqfil) hvis ingen subtypecalls
  # if("try-error" %in% class(RevReferences)){
  #   novcfcounter <- novcfcounter + 1
  #   next
  # }
  
  RevReferences <- RevReferences[,1]
  
  # Vælg om der kun skal bruges Most possible subtype (Fra Virstrain call) eller alle skal benyttes. Hvis alle benyttes kan det give et 
  # skævt billede af realiteten er, da hver fastqfil kan give flere resultater (nuværende 3, hvis ikke ændret) fra 1 SNV.
  # Hvis alle skal bruges, sættes 1:length(RevReferences ind istedet for 1 i følgende loop.
  # Måske mere relevant og eneste måde at fortsætte til nocalls script ligenu
  # er i stedet at tvinge en reference på alle fastqfiler. Sådan at eks. alle bliver tjekket imod K02718.1
  
  for(RRef in 1){

    ##TEST
    #RRef <- 1
    #Refname <- "K02718.1_revised"
    ######
    
    CallPrio <- which(RevReferences %in% RevReferences[RRef]) # Gør at det kan ses hvilken rang referencen har fra VirStrain
    counter <- counter + 1 # Tæller for at holde styr på hvor langt script er nået  
    Refname <- RevReferences[RRef]
    vcfname <- paste(Refname,".sort.dup.readGroupFix_filtered_FiltEx_headerfix.vcf", sep ="") # Tager fat i vcf med filteret varianter exluderet. # _filtered.filtEx_headerfix
    gffname <- paste(Refname,".gff3", sep ="")
    
    gfffile <- paste(GFFFolder, gffname, sep ="")
    Gene_anno <- makeTxDbFromGFF(gfffile, format = "gff3") # Laver TxDb objekt fra gff3 eller gtf fil. # , , circ_seqs = gfffile[1])
    
    Ref <- paste("/home/pato/Skrivebord/HPV16_projekt/References/IndexedRef/",Refname,"/", Refname, ".fasta", sep = "")
    faf <- open(FaFile(Ref))
    Folder <- paste("/home/pato/Skrivebord/HPV16_projekt/Results/",SuperRunName,"/", Runname,"/",Refname,"/", sep="")
    # Error check for om der er en filtreret vcf fil tilgængelig
    vcffile <- paste(Folder, vcfname, sep ="")
    
    c_vcf <- try(readVcf(vcffile)) # læser som collapsed vcf
    if("try-error" %in% class(c_vcf)){
      novcfcounter <- novcfcounter + 1
      next
    }
    
    # Error check som fortsætter til næste iteration i loopet hvis der ingen varianter i vcf fil er
    options(warn=2)
    no_elem <- try(predictCoding(c_vcf, Gene_anno, seqSource = faf))
    if("try-error" %in% class(no_elem)){
      noelement <- noelement + 1
      AnnoRdy <- data.frame(matrix(ncol = 1, nrow = 1))
      colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = "_")
      newdf <- cbind.fill(newdf, AnnoRdy)
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
      AnnoRdy[j,] <- paste(rdyNuc[j,],rdyAA[j,],Refname,sep = "%") 
    }
    
    colnames(AnnoRdy) <- paste(Fastqname,sep = "_")
    
    newdf <- cbind.fill(newdf, AnnoRdy)

    #newdf_nuc <- cbind.fill(newdf_nuc, rdyNuc)
    
    #newdf_nuc_coord <- cbind.fill(newdf_nuc_coord,Nuc$Nuc_start)
    
    # Conveniece for tracking progress
    print(paste(Fastqname,"is done...",counter,"of",lengthofList, sep = " "))
    
  }
  
  print("Job's done!")
  print(paste("Number of vcf files with no variants: ", noelement))
  print(paste("Number of vcf files corrupt or missing: ", novcfcounter))
  
}
  
# Find each unique instance of change:
Fredf <- table(newdf)
# Tester om der er nogen elementer
if(nrow(Fredf)==0) {
  Frequency <- data.frame(matrix(ncol=1,nrow=0))
  Fredf <- data.frame(matrix(ncol=1,nrow=0))
  df <- data.frame(matrix(ncol=1,nrow=0))
} else {
  Fredf <- as.data.frame(Fredf)
  # Splitting data table into separate columns
  Fredf <- separate(Fredf,col = 1, sep = "_", extra = "merge", into = c("Change", "GENEID&REF")) # Extra = merge sørger for at _ efter det første ikke bliver splittet
  Fredf <- separate(Fredf, col = "Change", sep = ",", extra = "merge", into = c("NucPos","NucChange&AA"))
  Fredf <- separate(Fredf, col = "NucChange&AA", sep = "%", extra = "merge", into = c("NucChange","AA"))
  Fredf <- separate(Fredf, col = "GENEID&REF", sep = "%", extra = "merge", into = c("GENEID","Reference"))
  Fredf$NucPos <- as.integer(Fredf$NucPos)
  Fredf <- Fredf[order(Fredf$Reference,Fredf$NucPos),]
}

# Laver en bed fil opsætning med nucpos for hver reference
Refs <- unique(Fredf$Reference)
for(cRef in Refs){
  cFredf <- Fredf[ Fredf$Reference == cRef, ]
  NucChangePos <- cFredf$NucPos
  df <- as.vector(as.matrix(NucChangePos))
  df <- data.frame(Reference = Refname, Pos = NucChangePos-1, PosEnd = NucChangePos)
  write.table(df, file = paste(SaveDir, "/",MultiFQfile,"_Nuc_change_coords_", cRef, ".bed", sep = ""), row.names = F,col.names = F, quote = F, sep = "\t")
}  

write.table(Fredf, file = paste(SaveDir, "/","ForNoCallScript_",MultiFQfile,"_Nuc_change_coords_", Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F, sep = "\t")

# Laver en pænere opsætning til tabel med individuelle
RdyDF <- newdf
for(c in 1:ncol(RdyDF)){
  for(r in 1:nrow(RdyDF)){
    
    RdyDF[r,c] <- str_replace_all(RdyDF[r,c], "%", " ")
    RdyDF[r,c] <- str_replace(RdyDF[r,c], "_", " ")
    RdyDF[r,c] <- str_remove(RdyDF[r,c], ",")
    RdyDF[r,c] <- paste("c.",RdyDF[r,c], sep ="")
    RdyDF[r,c] <- str_remove(RdyDF[r,c], "c.NA")
    
  }
}


write.table(RdyDF, file = paste(SaveDir, "/", "Annotations_",MultiFQfile, Refname, ".txt", sep = ""), row.names = F,col.names = T, quote = F, sep = "\t")



