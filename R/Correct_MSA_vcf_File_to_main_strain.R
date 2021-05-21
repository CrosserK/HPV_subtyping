
# THIS SCRIP:
# Correct_MSA_vcf_File_to_main_strain

# Find all positions of char in string and correct VCF positions
library(Biostrings)
library(tidyverse)

# Indsæt filer over ### og kør alt!
  MSA_afa <- "/home/pato/Skrivebord/HPV16_projekt/15_HPV_sublineages_aligned.afa" # MSA = Multiple sequence alignment
MSA_vcf <- "/home/pato/Skrivebord/HPV16_projekt/15_HPV_sublineages_aligned.vcf"   # MSA = Multiple sequence alignment
ChromosomeName <- "HQ644236.1" # Her kan navn ændres så filen kan læses i IGV. Det skal være tilsvarende referencen der loades i IGV
ID_header <- "##contig=<ID=1,length=7916>" # Her korrigeres til main strains længde
#################################

infile <- readDNAStringSet(MSA_afa)
# Dette vælger første sekvens i fasta: infile[1]
charpos <- gregexpr(pattern ='-',infile[1])
l <- charpos
gap_pos <- data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))
#gap_pos viser nu eksakte koordinater på "-"'er
# read two times the vcf file, first for the columns names, second for the data
vcf_header<-readLines(MSA_vcf)
vcf_data <- read.table(MSA_vcf, stringsAsFactors = FALSE)
tmp_vcf_data <- vcf_data[2]
names(tmp_vcf_data) <- "Pos"

#### test
startpos <- 1
endidx <- nrow(tmp_vcf_data)
for(j in gap_pos){
  highestpos <- j
  #print(highestpos)
  
  for(pos in tmp_vcf_data[,1]){
    #print(pos)
    if(pos >= highestpos){
      fromidx <- min(which(tmp_vcf_data[,1] >= highestpos))
      break()
    }
  }
  tmp_vcf_data[fromidx:endidx,1] <- tmp_vcf_data[fromidx:endidx,1] - 1
}
#####

vcf_data[,2] <- tmp_vcf_data[,1]
#view(vcf_header[1:4])
vcf_header[2] <- ID_header # Ændrer ifht. hvor mange gap_pos der er 
vcf_header <- vcf_header[1:4]
vcf_data[1] <- ChromosomeName # Korrigerer CHROM navn, så IGV kan læse den
write.table(vcf_header, file = "/home/pato/Skrivebord/HPV16_projekt/HPV_all_corr.vcf", sep = "\t", col.names = F, row.names = F, append = F, quote = F)
write.table(vcf_data, file = "/home/pato/Skrivebord/HPV16_projekt/HPV_all_corr.vcf", sep = "\t", col.names = F, row.names = F, append = T, quote = F)


# DER KOMMER TIL AT VÆRE NOGEN SOM HAR SAMME POSITION, DE SKAL OGSÅ KORRIGERES, 
# så de kommer på samme linje


