#!/usr/bin/env Rscript

# Predict amino acid changes from vcf, gff3 and fasta files, this version will output format for multiple fastqfiles
# in a convenient way

#install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("GenomicFeatures", force = TRUE)
#BiocManager::install("Repitools")
suppressMessages(library(GenomicFeatures))
suppressMessages(library(tidyverse))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(Repitools))
suppressMessages(library(dplyr))

# Henter variabler fra commandline argumenter
MainF <- as.character(commandArgs(TRUE)[1])
SaveDir <- as.character(commandArgs(TRUE)[2])
TopRunName <- as.character(commandArgs(TRUE)[3])
MultiFQfile <- as.character(commandArgs(TRUE)[4])
typeTier <- as.character(commandArgs(TRUE)[5])
LogFile <- paste(MainF,"/QC/Logs/",TopRunName,".txt", sep = "")


# MainF <- "/home/pato/Skrivebord/HPV_subtyping"
# TopRunName <- "Opsaet9_chip1_HPV16panel_0934_29112022"
# MultiFQfile <- paste("/home/pato/Skrivebord/HPV_subtyping/FASTQ/FASTQfiles_",TopRunName,".txt",sep="")
# SaveDir <- paste("/home/pato/Skrivebord/HPV_subtyping/Results/",TopRunName,sep="")
# typeTier <- "2"

# print(MainF)
# print(SaveDir)
# print(TopRunName)
# print(MultiFQfile)
# print(typeTier)

# Define function that enables columns of different length to be bound to each other, by filling shorter columns with "NA"
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
  rbind(x, matrix(,n-nrow(x), ncol(x))))) 
}

GFFFolder <- paste(MainF,"/References/GFFfiles/", sep="")

MultiFastqListFile <- MultiFQfile
MultiFastqList <- read.table(MultiFastqListFile)
MultiFastqList <- as.list(MultiFastqList)
MultiFastqList <- unlist(MultiFastqList)
splFolder <- paste(MainF,"/References/E4fixCoords/",sep = "")
newdf <- data.frame()
newdf_nuc <- data.frame()
newdf_nuc_coord <- data.frame()
lengthofList <- length(MultiFastqList)
counter <- 0
noelement <- 0
novcfcounter <- 0
for(Fastqname in MultiFastqList){
  #Fastqname <- MultiFastqList[1]
  #Fastqname <-  "/home/pato/Skrivebord/HPV_subtyping/FASTQ/Opsaet9_chip1_HPV16panel_1320_25112022/151.IonXpress_019.fastq" 
  Fastqname <- tools::file_path_sans_ext(basename(Fastqname)) # Getting fastqname without path and ext
  # Checking if there's split bamfiles because of multple HPV types in one fastq
  RevReferences <- try(read.table(paste(MainF,"/Results/",TopRunName,"/",Fastqname,"/TypeCalls/",Fastqname,"_T",typeTier,"_R1_SplitTo.txt", sep = "")))
  if("try-error" %in% class(RevReferences)){
    RevReferences <- try(read.table(paste(MainF,"/Results/",TopRunName,"/",Fastqname,"/TypeCalls/",Fastqname,"_T",typeTier,"_R1.txt", sep = "")))
    # If no reference (no type found) for fastq, go to next fastq
    if("try-error" %in% class(RevReferences)){
      next
    }
  }
  
  counter <- counter + 1
  
  # Vælg om der kun skal bruges Most possible subtype (Fra Virstrain call) eller alle skal benyttes. Hvis alle benyttes kan det give et 
  # skævt billede af realiteten er, da hver fastqfil kan give flere resultater (nuværende 3, hvis ikke ændret) fra 1 SNV.
  # Hvis alle skal bruges, sættes 1:length(RevReferences ind istedet for 1 i følgende loop.
  for(RRef in 1:1){
    # Enables checking which rank the reference has from  VirStrain, or given reference list
    # CallPrio <- which(RevReferences %in% RevReferences[RRef,]) 
    #RRef <- 1
    Refname <- RevReferences[RRef,]
    Folder <- paste(MainF,"/Results/",TopRunName,"/", Fastqname,"/",Refname,"/ResultFiles/", sep="")
    # Checking if there is any vcf in the ref folder (if not, probably nothing found in bamsplit)
    
    vcfname <- paste(Fastqname,"_", Refname,"_filt.vcf", sep = "")
    gffname <- paste(Refname,".gff3", sep ="")
    gfffile <- paste(GFFFolder, gffname, sep ="")
    print(paste("Using GFF",gffname,"for",Fastqname,sep=" "))
    Gene_anno <- makeTxDbFromGFF(gfffile, format = "gff3") # Laver TxDb objekt fra gff3 eller gtf fil. # , , circ_seqs = gfffile[1])
    
    Ref <- paste(MainF,"/References/IndexedRef/",Refname,"/", Refname, ".fasta", sep = "")
    faf <- open(FaFile(Ref))
    
    # Error check for om der er en filtreret vcf fil tilgængelig
    vcffile <- paste(Folder, vcfname, sep ="")
    c_vcf <- try(readVcf(vcffile)) # læser som collapsed vcf
    # Check that there is variants in the vcf file, else go to next reference
    
    # Checker om c_vcf har info som den skal
    if(length(dim(c_vcf))>1){
      if(dim(c_vcf)[1] < 1){
        noelement <- noelement + 1
        next
      }
    } else {
      noelement <- noelement + 1
      next
    }
    
    if("try-error" %in% class(c_vcf)){
      novcfcounter <- novcfcounter + 1
      AnnoRdy <- data.frame(matrix(ncol = 1, nrow = 1))
      colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
      AnnoRdy[1] <- "NoVcfFile"
      newdf <- cbind.fill(newdf, AnnoRdy)
      next
    }
    
    # Error check som fortsætter til næste iteration i loopet hvis der ingen varianter i vcf fil er. Der er en anden error message, 
    # fra predictCoding, som kan fortælle at varianter er out-of-bounds, men denne virker til at være bugged, når gff fil har flere exons.
    # Derfor ignoreres den og der ledes kun efter fejl med string "none of [..]", som nedenunder
    vcf_anno <- try(predictCoding(c_vcf, Gene_anno, seqSource = faf)) # Læser som collapsed vcf
    
    if("try-error" %in% class(vcf_anno)){
      noelement <- noelement + 1
      AnnoRdy <- data.frame(matrix(ncol = 1, nrow = 1))
      colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
      AnnoRdy[1] <- "NO ELEMENTS IN VCF FILE"
      newdf <- cbind.fill(newdf, AnnoRdy)
      next
    }
    
    Nuc_start <- vcf_anno@ranges@start # Collecting nuc position
    
    vcf_anno_df_2 <- vcf_anno@elementMetadata # Collecting elementdata
    vcf_anno_df <- cbind(Nuc_start, vcf_anno_df_2)
    anno_df <- as.data.frame(vcf_anno_df)
    
    # Insert missing amino acid, when frameshift
    for(i in 1:nrow(anno_df)) {
      if(anno_df$REFAA[i] == ""){
        codon <- DNAString(anno_df$REFCODON[i])
        codonTra <- suppressWarnings(Biostrings::translate(codon))
        anno_df$REFAA[i] <- as.character(codonTra) 
      } 
    }
    

    # BEGYNDER E4 SPLICE GEN KORREKTION:
    # Henter faktiske koordinater fra spliceinfo gff fil
    splname <- paste(Refname,".gff3", sep ="")
    splfile <- paste(splFolder, splname, sep ="")
    
    # Finder splice koordinater
    print(paste("Using splice file",splfile))
    txdb <- try(makeTxDbFromGFF(splfile, format = "gff3"))
    if("try-error" %in% class(txdb)){
      AnnoRdy <- data.frame(matrix(ncol = 1, nrow = 1))
      colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
      AnnoRdy[1] <- paste("E4FixGFF3DoesNotExistOrIsCorrupt;",Refname,sep="")
      newdf <- cbind.fill(newdf, AnnoRdy)
      next
    }
    txseq <- extractTranscriptSeqs(faf,cdsBy(txdb, by="tx", use.names=TRUE))
    translated <- suppressWarnings(Biostrings::translate(txseq))
    gffDF <- as.data.frame(VariantAnnotation::select(txdb, column=columns(txdb), keys="E4_splice", keytype=("GENEID")))  # Fanger gff fil info for at hente koordinater på splicetranskript
    
    # For fat på alt nucs i liste
    altNuc <- anno_df$ALT
    altNuc = lapply(altNuc, as.character)
    altNuc = sapply(altNuc, paste, collapse = ",")
    # Oversætter nuc koord til splice koord: 
    CoordsFirstPart <- c(gffDF$CDSSTART[1],gffDF$CDSEND[1])  # Få fat på faktisk codon
    # Der skal +1 fordi ellers tælles alle nuc ikke med
    FirstSplceLen <- CoordsFirstPart[2] - CoordsFirstPart[1] + 1
    CoordsSecondPart <- c(gffDF$CDSSTART[2],gffDF$CDSEND[2])
    # Minus 1 for ellers tælles en nucleotid for meget i intron længde
    
    spliceseq <- as.character(txseq$E4_spliceRNA)
    splceAA <- as.character(translated$E4_spliceRNA)
    
    # For hver variant, tjek for ukorrigeret splicegen kald
    for(annoRow in 1:nrow(anno_df)){
      
      #annoRow <- 14
      # Check if row exists
      rowExistsCheck <- try(anno_df$GENEID[annoRow])
      
      # Check if row exists. If NA, set the variant as not in gene.
      if(length(rowExistsCheck) == 0 || is.na(rowExistsCheck)){
        if(exists("AnnoRdy")){
          colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
          AnnoRdy[1] <- "ErrorCannotHandleGff"
          #newdf <- cbind.fill(newdf, AnnoRdy)
          next
        } else {
          AnnoRdy <- data.frame(matrix(ncol = 1, nrow = 1))
          colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
          AnnoRdy[1] <- "ErrorCannotHandleGff"
          #newdf <- cbind.fill(newdf, AnnoRdy)
          next
        }
        
      } else {
        # If it is E4_splice
        if(anno_df$GENEID[annoRow] == "E4_splice"){
        # Konverterer DNAstringsets til characters
        altAsChar <- sapply(anno_df$ALT[annoRow], function(x) as.character(x[[1]]))
        refAsChar <- sapply(anno_df$REF[annoRow], function(x) as.character(x[[1]]))
        
        # Hvis længden af variant = 1:
        if(nchar(refAsChar) == 1 && nchar(altAsChar) == 1){
          print("E4 SNV")
          
          # Finder koordinat i splice gen
          if(anno_df$Nuc_start[annoRow] <= CoordsFirstPart[2]){
            # Henter position, derfor +1 (eks hvis pos 892, svarer til pos 1, så minuses 892 med 892 og +1)
            splicecoord <- anno_df$Nuc_start[annoRow] - CoordsFirstPart[1] + 1  
          } else if(anno_df$Nuc_start[annoRow] <= CoordsSecondPart[2]){
            splicecoord <- anno_df$Nuc_start[annoRow] + FirstSplceLen - (CoordsSecondPart[1]-1)
          }
          
          
          # Finder nukleotides position i codon
          invPosInCodon <- splicecoord%%3 
          # Finder codon:
          if(invPosInCodon == 0){
            codon <- substring(spliceseq,splicecoord-2,splicecoord)
            PosInCodon <- 2
          } else if(invPosInCodon == 1){
            codon <- substring(spliceseq,splicecoord,splicecoord+2)
            PosInCodon <- 0
          } else if(invPosInCodon == 2){
            codon <- substring(spliceseq,splicecoord-1,splicecoord+1)
            PosInCodon <- 1
          }
          
          ThirdNucPos <- 2 - PosInCodon
          # Retter kald til faktisk splicet gen
          AAPos <- (splicecoord + ThirdNucPos)/3  # Finding amino acid
          actualAA <- substr(splceAA,AAPos,AAPos)
          
          # Bruger integreret AA converter ved at give fuld sekvens
          currentAltNuc <- altNuc[annoRow]
          
          seqForAAconv <- spliceseq
          substr(seqForAAconv, splicecoord, splicecoord) <- currentAltNuc
          
          # Finder ny alt codon
          if(invPosInCodon == 0){
            altCodon <- substring(seqForAAconv,splicecoord-2,splicecoord)
          } else if(invPosInCodon == 1){
            altCodon <- substring(seqForAAconv,splicecoord,splicecoord+2)
          } else if(invPosInCodon == 2){
            altCodon <- substring(seqForAAconv,splicecoord-1,splicecoord+1)
          }
          
          # Translaterer
          seqForAAconv <- DNAStringSet(seqForAAconv)
          actualTrans <- suppressWarnings(Biostrings::translate(seqForAAconv))
          
          # Laver tilbage til string:
          actualTransChar <- as.character(actualTrans)
          actualALTAA <- substring(actualTransChar,AAPos, AAPos)
          
          # Ellers hvis der er fundet en deletion
        } else if(nchar(refAsChar) > 1 && nchar(altAsChar) == 1){
          print("E4 del")
          
          # Tjekker for frameshift:
          if((nchar(refAsChar)-1)%%3 == 0){
            
            # Ingen frameshift, finder aminosyre ændring
            # Finder koordinat i splice gen
            if(anno_df$Nuc_start[annoRow] <= CoordsFirstPart[2]){
              # Henter position, derfor +1 (eks hvis pos 892, svarer til pos 1, så minuses 892 med 892 og +1)
              splicecoord <- anno_df$Nuc_start[annoRow] - CoordsFirstPart[1] + 1  
            } else if(anno_df$Nuc_start[annoRow] <= CoordsSecondPart[2]){
              splicecoord <- anno_df$Nuc_start[annoRow] + FirstSplceLen - (CoordsSecondPart[1]-1)
            }
            
            extraNCodons <- (nchar(refAsChar)-1)/3
            
            # Finder nukleotides position i codon og codon
            invPosInCodon <- splicecoord%%3 
            # Finder codon:
            
            if(invPosInCodon == 0){
              firstCodon <- substring(spliceseq,splicecoord-2,splicecoord)
              PosInCodon <- 2
              codon <- paste(firstCodon, substr(spliceseq,splicecoord+1,splicecoord+1+extraNCodons*3-1), sep="")
            } else if(invPosInCodon == 1){
              firstCodon <- substring(spliceseq,splicecoord,splicecoord+2)
              PosInCodon <- 0
              codon <- paste(firstCodon, substr(spliceseq,splicecoord+3,splicecoord+3+extraNCodons*3-1), sep="")
            } else if(invPosInCodon == 2){
              firstCodon <- substring(spliceseq,splicecoord-1,splicecoord+1)
              PosInCodon <- 1
              codon <- paste(firstCodon, substr(spliceseq,splicecoord+2,splicecoord+2+extraNCodons*3-1), sep="")
            }
            
            ThirdNucPos <- 2 - PosInCodon
            # Retter kald til faktisk splicet gen
            # Finder AA:
            AAPos <- (splicecoord + ThirdNucPos)/3 
            actualAA <- paste("c(",AAPos,":",AAPos+extraNCodons,")", substr(splceAA,AAPos,AAPos + extraNCodons),sep="")
            
            # Bruger integreret aa converter ved at give fuld sekvens
            currentAltNuc <- altNuc[annoRow]
            
            seqForAAconv <- spliceseq
            substr(seqForAAconv, splicecoord, splicecoord) <- currentAltNuc
            
            
            # Finder ny alt codon
            if(invPosInCodon == 0){
              altCodon <- substring(seqForAAconv,splicecoord-2,splicecoord)
            } else if(invPosInCodon == 1){
              altCodon <- substring(seqForAAconv,splicecoord,splicecoord+2)
            } else if(invPosInCodon == 2){
              altCodon <- substring(seqForAAconv,splicecoord-1,splicecoord+1)
            }
            
            # Translaterer
            seqForAAconv <- DNAStringSet(seqForAAconv)
            actualTrans <- suppressWarnings(Biostrings::translate(seqForAAconv))
            
            # Laver tilbage til string:
            actualTransChar <- as.character(actualTrans)
            actualALTAA <- substring(actualTransChar,AAPos, AAPos) 
            
          } else {
            
            # Frameshift
            actualALTAA <- "*"
            
            # Ingen frameshift, finder aminosyre ændring
            # Finder koordinat i splice gen
            if(anno_df$Nuc_start[annoRow] <= CoordsFirstPart[2]){
              # Henter position, derfor +1 (eks hvis pos 892, svarer til pos 1, så minuses 892 med 892 og +1)
              splicecoord <- anno_df$Nuc_start[annoRow] - CoordsFirstPart[1] + 1  
            } else if(anno_df$Nuc_start[annoRow] <= CoordsSecondPart[2]){
              splicecoord <- anno_df$Nuc_start[annoRow] + FirstSplceLen - (CoordsSecondPart[1]-1)
            }
            
            extraNCodons <- (nchar(refAsChar)-1)/3
            
            # Finder nukleotides position i codon og codon
            invPosInCodon <- splicecoord%%3 
            # Finder codon:
            
            if(invPosInCodon == 0){
              firstCodon <- substring(spliceseq,splicecoord-2,splicecoord)
              PosInCodon <- 2
              codon <- paste(firstCodon, substr(spliceseq,splicecoord+1,splicecoord+1+extraNCodons*3-1), sep="")
            } else if(invPosInCodon == 1){
              firstCodon <- substring(spliceseq,splicecoord,splicecoord+2)
              PosInCodon <- 0
              codon <- paste(firstCodon, substr(spliceseq,splicecoord+3,splicecoord+3+extraNCodons*3-1), sep="")
            } else if(invPosInCodon == 2){
              firstCodon <- substring(spliceseq,splicecoord-1,splicecoord+1)
              PosInCodon <- 1
              codon <- paste(firstCodon, substr(spliceseq,splicecoord+2,splicecoord+2+extraNCodons*3-1), sep="")
            }
            
            ThirdNucPos <- 2 - PosInCodon
            # Retter kald til faktisk splicet gen
            # Finder AA:
            AAPos <- (splicecoord + ThirdNucPos)/3 
            actualAA <- paste("c(",AAPos,":",AAPos+extraNCodons,")", substr(splceAA,AAPos,AAPos + extraNCodons),sep="")
            
            # Bruger integreret aa converter ved at give fuld sekvens
            currentAltNuc <- altNuc[annoRow]
            
            seqForAAconv <- spliceseq
            substr(seqForAAconv, splicecoord, splicecoord) <- currentAltNuc
            
            
            # Finder ny alt codon
            if(invPosInCodon == 0){
              altCodon <- substring(seqForAAconv,splicecoord-2,splicecoord)
            } else if(invPosInCodon == 1){
              altCodon <- substring(seqForAAconv,splicecoord,splicecoord+2)
            } else if(invPosInCodon == 2){
              altCodon <- substring(seqForAAconv,splicecoord-1,splicecoord+1)
            }
            
          }
          
          # Ellers hvis der er fundet en insertion
        } else if(nchar(refAsChar) == 1 && nchar(altAsChar) > 1){
          print("E4 ins")
          
          # Tjekker for frameshift:
          if((nchar(refAsChar)-1)%%3 == 0){
            
            # Ingen frameshift, finder aminosyre ændring
            # Finder koordinat i splice gen
            if(anno_df$Nuc_start[annoRow] <= CoordsFirstPart[2]){
              # Henter position, derfor +1 (eks hvis pos 892, svarer til pos 1, så minuses 892 med 892 og +1)
              splicecoord <- anno_df$Nuc_start[annoRow] - CoordsFirstPart[1] + 1  
            } else if(anno_df$Nuc_start[annoRow] <= CoordsSecondPart[2]){
              splicecoord <- anno_df$Nuc_start[annoRow] + FirstSplceLen - (CoordsSecondPart[1]-1)
            }
            
            
            extraNCodons <- (nchar(altAsChar)-1)/3
            
            # Finder nukleotides position i codon og finder codon
            invPosInCodon <- splicecoord%%3 
            # Finder codon:
            
            if(invPosInCodon == 0){
              codon <- substring(spliceseq,splicecoord-2,splicecoord)
              PosInCodon <- 2
            } else if(invPosInCodon == 1){
              codon <- substring(spliceseq,splicecoord,splicecoord+2)
              PosInCodon <- 0
            } else if(invPosInCodon == 2){
              codon <- substring(spliceseq,splicecoord-1,splicecoord+1)
              PosInCodon <- 1
            }
            
            ThirdNucPos <- 2 - PosInCodon
            # Retter kald til faktisk splicet gen
            # Finder AA:
            AAPos <- (splicecoord + ThirdNucPos)/3 
            actualAA <- substr(splceAA,AAPos,AAPos)
            
            # Bruger integreret aa converter ved at give fuld sekvens
            currentAltNuc <- altAsChar
            seqForAAconv <- spliceseq
            # Sætter ny splice sekvens sammen med insertion
            seqForAAconv <- paste(substr(seqForAAconv, 1, splicecoord),substr(currentAltNuc,2,nchar(currentAltNuc)),substr(seqForAAconv, splicecoord+1, nchar(seqForAAconv)),sep="") 
            
            # Finder ny alt codon
            if(invPosInCodon == 0){
              firstCodon <- substring(seqForAAconv,splicecoord-2,splicecoord)
              # Udvider codon med insertion
              altCodon <- paste(firstCodon,substring(seqForAAconv,splicecoord+1,splicecoord+1+extraNCodons*3-1),sep="")
            } else if(invPosInCodon == 1){
              firstCodon <- substring(seqForAAconv,splicecoord,splicecoord+2)
              # Udvider codon med insertion
              altCodon <- paste(firstCodon,substring(seqForAAconv,splicecoord+3,splicecoord+3+extraNCodons*3-1),sep="")
            } else if(invPosInCodon == 2){
              firstCodon <- substring(seqForAAconv,splicecoord-1,splicecoord+1)
              # Udvider codon med insertion
              altCodon <- paste(firstCodon,substring(seqForAAconv,splicecoord+2,splicecoord+2+extraNCodons*3-1),sep="")
            }
            
            # Translaterer
            seqForAAconv <- DNAStringSet(seqForAAconv)
            actualTrans <- suppressWarnings(Biostrings::translate(seqForAAconv))
            
            # Laver tilbage til string:
            actualTransChar <- as.character(actualTrans)
            actualALTAA <- substring(actualTransChar,AAPos, AAPos+extraNCodons) 
            
            
            } else {
              
              # Frameshift
              actualALTAA <- "*"
              
              # Ingen frameshift, finder aminosyre ændring
              # Finder koordinat i splice gen
              if(anno_df$Nuc_start[annoRow] <= CoordsFirstPart[2]){
                # Henter position, derfor +1 (eks hvis pos 892, svarer til pos 1, så minuses 892 med 892 og +1)
                splicecoord <- anno_df$Nuc_start[annoRow] - CoordsFirstPart[1] + 1  
              } else if(anno_df$Nuc_start[annoRow] <= CoordsSecondPart[2]){
                splicecoord <- anno_df$Nuc_start[annoRow] + FirstSplceLen - (CoordsSecondPart[1]-1)
              }
              
              
              extraNCodons <- (nchar(altAsChar)-1)/3
              
              # Finder nukleotides position i codon og codon
              invPosInCodon <- splicecoord%%3 
              # Finder codon:
              
              if(invPosInCodon == 0){
                codon <- substring(spliceseq,splicecoord-2,splicecoord)
                PosInCodon <- 2
              } else if(invPosInCodon == 1){
                codon <- substring(spliceseq,splicecoord,splicecoord+2)
                PosInCodon <- 0
              } else if(invPosInCodon == 2){
                codon <- substring(spliceseq,splicecoord-1,splicecoord+1)
                PosInCodon <- 1
              }
              
              ThirdNucPos <- 2 - PosInCodon
              # Retter kald til faktisk splicet gen
              # Finder AA:
              AAPos <- (splicecoord + ThirdNucPos)/3 
              actualAA <- substr(splceAA,AAPos,AAPos)
              
              # Bruger integreret aa converter ved at give fuld sekvens
              currentAltNuc <- altAsChar
              seqForAAconv <- spliceseq
              # Sætter ny splice sekvens sammen med insertion
              seqForAAconv <- paste(substr(seqForAAconv, 1, splicecoord),substr(currentAltNuc,2,nchar(currentAltNuc)),substr(seqForAAconv, splicecoord+1, nchar(seqForAAconv)),sep="") 
              
              # Finder ny alt codon
              if(invPosInCodon == 0){
                firstCodon <- substring(seqForAAconv,splicecoord-2,splicecoord)
                # Udvider codon med insertion
                altCodon <- paste(firstCodon,substring(seqForAAconv,splicecoord+1,splicecoord+1+extraNCodons*3-1),sep="")
              } else if(invPosInCodon == 1){
                firstCodon <- substring(seqForAAconv,splicecoord,splicecoord+2)
                # Udvider codon med insertion
                altCodon <- paste(firstCodon,substring(seqForAAconv,splicecoord+3,splicecoord+3+extraNCodons*3-1),sep="")
              } else if(invPosInCodon == 2){
                firstCodon <- substring(seqForAAconv,splicecoord-1,splicecoord+1)
                # Udvider codon med insertion
                altCodon <- paste(firstCodon,substring(seqForAAconv,splicecoord+2,splicecoord+2+extraNCodons*3-1),sep="")
              }
                
              
            }
            
        }
        # Retter tabel med splice informationer
        anno_df$PROTEINLOC[annoRow] <- AAPos
        anno_df$QUERYID[annoRow] <- "x"
        anno_df$TXID[annoRow] <- "x"
        anno_df$CDSID[annoRow] <- "x"
        #anno_df$CONSEQUENCE[annoRow] <- as.factor("x")
        anno_df$REFCODON[annoRow] <- codon
        anno_df$VARCODON[annoRow] <- altCodon
        anno_df$REFAA[annoRow] <- actualAA
        anno_df$VARAA[annoRow] <- actualALTAA
        }
      }
      
      # Formatting protein changes
      # Checking that variants are inside genes of gff file
      AAchange <- try(anno_df[c("REFAA","PROTEINLOC", "VARAA","GENEID")]) # "REFAA","PROTEIONLOC","VARAA")
      if("try-error" %in% class(AAchange)){
        print(paste("No variants in genes of GFF file",Refname))
        rdyAA <- as.data.frame(paste("VarOutsideGeneRegionOfGFF;",Refname, sep = ""), col.names = aa)
        rdyNuc <- as.data.frame(paste("VarOutsideRegionOfGFF;",Refname, sep = ""))
        colnames(rdyNuc) <- "Mutation"

        #Merger AA og nuc info til 1 kolonne
        AnnoRdy <- data.frame(matrix(ncol = 1, nrow = nrow(rdyAA)))
        for(j in 1:nrow(rdyAA)){
          AnnoRdy[j,] <- paste("VarOutsideRegionOfGFF;",Refname, sep = "") 
        }
        colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
        newdf <- cbind.fill(newdf, AnnoRdy)
        next
        
        # If there are variants inside gene regions of gff file
      } else {
        
        # Changing empty varaa and refaa to frameshift
        AAchange['VARAA'][AAchange['VARAA'] == ""] <- "Frameshift"
        
        AA <- cbind(p = "p.", AAchange)
        rdyAA <- as.data.frame(paste(AA$p, AA$REFAA, AA$PROTEINLOC, AA$VARAA,":", AA$GENEID, sep = ""), col.names = aa)
        
        # Formatting nucleotide changes
        Nuchange <- vcf_anno_df[c("Nuc_start","REF", "ALT")] # "REFAA","PROTEIONLOC","VARAA"
        Nuc <- cbind(Nuchange, a = ">")
        #Nuc <- cbind(c = "c.", Nuc)
        Nuc <- cbind(sepera = ",", Nuc)
        rdyNuc <- as.data.frame(paste(Nuc$Nuc_start,Nuc$sepera, Nuc$REF, Nuc$a, Nuc$ALT@unlistData, sep = ""))
        colnames(rdyNuc) <- "Mutation"
        #Merger AA og nuc info til 1 kolonne
        AnnoRdy <- data.frame(matrix(ncol = 1, nrow = nrow(rdyAA)))
        for(j in 1:nrow(rdyAA)){
          AnnoRdy[j,] <- paste(rdyNuc[j,],rdyAA[j,],Refname,sep = "%") 
        }
      }
    }
        
 
    # Convenience for tracking progress
    #print(paste(Fastqname,"is done...",counter,"of",lengthofList, sep = " "))
    print(paste("Number of vcf files with no variants: ", noelement))
    print(paste("Number of vcf files corrupt or missing (possibly no bam from bamsplit): ", novcfcounter))  
    colnames(AnnoRdy) <- paste(Fastqname,Refname,sep = ":")
    newdf <- cbind.fill(newdf, AnnoRdy)
    #newdf_nuc <- cbind.fill(newdf_nuc, rdyNuc)
    #newdf_nuc_coord <- cbind.fill(newdf_nuc_coord,Nuc$Nuc_start)
  }
}

# For log
text1 <- paste("Number of vcf files with no variants: ", noelement)
text2 <- paste("Number of vcf files corrupt or missing (possibly no bam from bamsplit): ", novcfcounter)
write.table(text1, file = LogFile, append = TRUE, quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
write.table(text2, file = LogFile, append = TRUE, quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)


# Find each unique instance of change:
# Tester om der er nogen elementer
if(nrow(newdf)==0) {
  Frequency <- data.frame(matrix(ncol=1,nrow=0))
  Fredf <- data.frame(matrix(ncol=1,nrow=0))
  df <- data.frame(matrix(ncol=1,nrow=0))
} else {
  Fredf <- table(newdf)
  Fredf <- as.data.frame(Fredf)
  # Splitting data table into separate columns
  Fredf <- separate(Fredf, col = 1, sep = ":", extra = "merge", into = c("Change", "GENEID&REF"), fill = "right") # Extra = merge sørger for at _ efter det første ikke bliver splittet
  Fredf <- separate(Fredf, col = "Change", sep = ",", extra = "merge", into = c("NucPos","NucChange&AA"), fill = "right")
  Fredf <- separate(Fredf, col = "NucChange&AA", sep = "%", extra = "merge", into = c("NucChange","AAChange"), fill = "right")
  Fredf <- separate(Fredf, col = "GENEID&REF", sep = "%", extra = "merge", into = c("GENEID","Reference"), fill = "right")
  # Moving error message from nucPos to GENEID
  for(fredfRow in 1:nrow(Fredf)){
    if(is.na(Fredf$GENEID[fredfRow])){
      Fredf$GENEID[fredfRow] <- Fredf$NucPos[fredfRow]
    }
  }
  Fredf$NucPos <- as.integer(Fredf$NucPos)
  Fredf <- Fredf[order(Fredf$Reference,Fredf$NucPos),]
}



# Laver en bed fil opsætning med nucpos for hver reference
Refs <- unique(Fredf$Reference[!is.na(Fredf$Reference)])
Refs

for(cRef in Refs){
  cFredf <- Fredf[ Fredf$Reference == cRef & !is.na(Fredf$Reference) , ]
  NucChangePos <- cFredf$NucPos[!is.na(cFredf$NucPos)]
  NucChangePos <- as.numeric(NucChangePos)
  df <- as.vector(as.matrix(NucChangePos))
  df <- data.frame(Reference = cRef, Pos = NucChangePos-1, PosEnd = NucChangePos)
  write.table(df, file = paste(SaveDir, "/",TopRunName,"_Nuc_change_coords_", cRef, ".bed", sep = ""), row.names = F,col.names = F, quote = F, sep = "\t")
}

write.table(Fredf, file = paste(SaveDir, "/","ForSummaryScript_",TopRunName,"_Nuc_change_coords.txt", sep = ""), row.names = F,col.names = T, quote = F, sep = "\t")

# Laver en pænere opsætning til tabel med individuelle
RdyDF <- newdf

if(length(RdyDF > 0)){
  for(c in 1:ncol(RdyDF)){
    for(r in 1:nrow(RdyDF)){
      if(!is.na(RdyDF[r,c])){
        if(RdyDF[r,c]!="NoVcfFile"){
            RdyDF[r,c] <- str_replace_all(RdyDF[r,c], "%", " ")
            RdyDF[r,c] <- str_replace(RdyDF[r,c], ":", " ")
            RdyDF[r,c] <- str_remove(RdyDF[r,c], ",")
            RdyDF[r,c] <- paste("c.",RdyDF[r,c], sep ="")
            RdyDF[r,c] <- str_remove(RdyDF[r,c], "c.NA")
        } 
      }
    }
  }
  RdyDF <- as.data.frame(RdyDF)
  
  # Fjerner tomme kolonner
  # RdyDF <- RdyDF[, colSums(RdyDF != "") != 0, drop = F] # drop = F sørger for at dataframe ikke taber kolonnenavne ved subsetting
  
  write.table(RdyDF, file = paste(SaveDir, "/", "AnnotationIndividualFiles_",TopRunName,".txt", sep = ""), na = "" , row.names = F,col.names = T, quote = F, sep = "\t")
  
} else {
  print("No variants found. Low coverage?")
}

