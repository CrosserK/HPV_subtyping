
# THIS SCRIP:
# Correct_MSA_vcf_File_to_main_strain

# Find all positions of gaps ('-') in string and correct VCF positions. Does not output gap positions, because SNP-sites does not.
library(Biostrings)
library(tidyverse)
library(sjmisc)


# Indsæt filer over ### og kør alt!
MainF <- as.character(commandArgs(TRUE)[1])
MSA_afa_navn <- as.character(commandArgs(TRUE)[2])
vcfnavn <- as.character(commandArgs(TRUE)[3]) # Inkl extenstion 
ChromosomeName <- as.character(commandArgs(TRUE)[4]) # Her kan navn ændres så filen kan læses i IGV. Det skal være tilsvarende referencen der loades i IGV
ID_header <- "##contig=<ID=1,length=na>" # Her korrigeres til main strains længde
#################################

MainF <- "/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types"
MSA_afa_navn <- "HPV97.mafft"
vcfnavn <- "HPV97.vcf" # Inkl extenstion
ChromosomeName <- "DQ080080.1" # Her kan navn ændres så filen kan læses i IGV. Det skal være tilsvarende referencen der loades i IGV
ID_header <- "##contig=<ID=1,length=na>"

MainFvcf <- paste(MainF, "/SNPs/",sep="")
MainFmafft <- paste(MainF, "/MSA/",sep="")

MSA_afa <- paste(MainFmafft,MSA_afa_navn,sep="") # MSA = Multiple sequence alignment
MSA_vcf <- paste(MainFvcf,vcfnavn, sep="")   # MSA = Multiple sequence alignment
infile <- readDNAStringSet(MSA_afa)
# Dette vælger første sekvens i fasta: infile[1]
charpos <- gregexpr(pattern ='-',infile[1]) # Finder alle gaps i øverste sekvens fra MSA
l <- charpos
vcf_header<-readLines(MSA_vcf)
vcf_header[2] <- ID_header # Ændrer ifht. hvor mange gap_pos der er 
vcf_header <- vcf_header[1:4]
vcf_data <- read.table(MSA_vcf, stringsAsFactors = FALSE)
if(l[[1]][1] > 0){ # Håndterer cases hvor der ingen gaps er
  gap_pos <- data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))
  ####TES
  #gap_pos[4] <- 3952
  #gap_pos viser nu eksakte koordinater på "-"'er
  # read two times the vcf file, first for the columns names, second for the data

  tmp_vcf_data <- vcf_data[2]
  names(tmp_vcf_data) <- "Pos"
  #tmp_vcf_data <- as.numeric(unlist(tmp_vcf_data))
  #tmp_vcf_data[4,1] <- 7578
  
  ####
  startpos <- 1
  endidx <- nrow(tmp_vcf_data)
  for(j in gap_pos){
    highestpos <- j
    #print(highestpos)
    
    for(pos in tmp_vcf_data[,1]){
      #print(pos)
      if(pos >= highestpos){
        fromidx <- min(which(tmp_vcf_data[,1] >= highestpos))
        tmp_vcf_data[fromidx:endidx,1] <- tmp_vcf_data[fromidx:endidx,1] - 1 # Trækker alle positioner efter gap sammen ved at minus dem med 1
        # Hopper ud af loop for denne gap pos og går til næste gap pos
        break()
      }
    }
  }
  
  #####
  vcf_data[,2] <- tmp_vcf_data[,1]
  #view(vcf_data)
  #view(vcf_data[,10:24])
  #rowSums(vcf_data[,10:24], x == "1")
  #tæl <- as.data.frame(vcf_data[,10:24])
  #row_count(tæl, count == "1")
  # #Indsætter nu gap positioner
  # Tager vilkårlig række til kopiering
  # vcf_data[1,]
  # for(r in gap_pos){
  #   
  # }
  # vcf_data <- rbind(vcf_data[1:r,],newrow,existingDF[-(1:r),])
  
  #str(vcf_data)
  #view(vcf_header[1:4])
  
  vcf_data[1] <- gsub("\\..*","",paste(substr(vcfnavn, 1, nchar(vcfnavn)-4),"_",ChromosomeName, sep="")) # Korrigerer CHROM navn, så IGV kan læse den. Korrigerer navn.
  write.table(vcf_header, file = paste(MainF,"/SNPs_corr/","corr_",vcfnavn,sep=""), sep = "\t", col.names = F, row.names = F, append = F, quote = F)
  write.table(vcf_data, file = paste(MainF,"/SNPs_corr/","corr_",vcfnavn,sep=""), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
  
  
  # DER KOMMER MULIGIVS TIL AT VÆRE NOGEN SOM HAR SAMME POSITION, DE SKAL OGSÅ KORRIGERES, 
  # så de kommer på samme linje
  
} else { # Skriver eksakt kopi ud, hvis der ingen gaps var i multiple seq alignment
  vcf_data[1] <- gsub("\\..*","",paste(substr(vcfnavn, 1, nchar(vcfnavn)-4),"_",ChromosomeName, sep=""))
  write.table(vcf_header, file = paste(MainF,"/SNPs_corr/","corr_",vcfnavn,sep=""), sep = "\t", col.names = F, row.names = F, append = F, quote = F)
  write.table(vcf_data, file = paste(MainF,"/SNPs_corr/","corr_",vcfnavn,sep=""), sep = "\t", col.names = F, row.names = F, append = T, quote = F)
}
# DER KOMMER MULIGIVS TIL AT VÆRE NOGEN SOM HAR SAMME POSITION, DE SKAL OGSÅ KORRIGERES, 
# så de kommer på samme linje


