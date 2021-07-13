library(tidyverse)
  
# Define folders
SuperRunName <- "Megakørsel1"
MainF <- "/home/pato/Skrivebord/HPV16_projekt"
##########

SigmaF <- paste(MainF,"/Sigma_run", sep = "")

# Get files. Finder Mismatch count filer i mapper navngivet på denne måde Pt_*_run 
Folderlist <- list.files(path = SigmaF, pattern = "^Pt_(.*)_run$", include.dirs = T)

# Laver start tabel
#data.frame <- 
  
Fastq1 <- paste(SigmaF,"/", Folder[1], "/MismatchCounts.txt", sep = "")
RefTable1 <- read.table(CurrentFastq)
SuperList <- data.frame(matrix(NA, nrow = nrow(RefTable1), ncol = 1))
colnames(SuperList) <- c("Ref")
SuperList[1] <- RefTable1[1:nrow(RefTable1),1]

for(Folder in Folderlist){
  CurrentFastq <- paste(SigmaF,"/", Folder, "/MismatchCounts.txt", sep = "")
  Currenttable <- read.table(CurrentFastq)
  colnames(Currenttable) <- c("Ref", Folder)
  SuperList <- merge(SuperList, Currenttable, by = "Ref", all = T)
}

row.names(SuperList) <- SuperList[1:nrow(SuperList),1]
SuperList <- SuperList[-1]

# Skalerer mismatches til højeste MM fundet i fastq

NormSuperList <- SuperList
for(n in 1:ncol(SuperList)){
  maxMM <- max(SuperList[n])
  NormSuperList[n] <- SuperList[n] / maxMM
}
view(NormSuperList)

# Fjerner kolonner som kun består af NA
NormSuperList <- NormSuperList[ , colSums(is.na(NormSuperList)) < nrow(NormSuperList)]

# Finder hyppigst forekommende subtype:
Ref_occur_tracker <- data.frame(matrix(0, nrow = nrow(NormSuperList), ncol = 1)) 
row.names(Ref_occur_tracker) <- row.names(NormSuperList)
colnames(Ref_occur_tracker) <- "Occurences"
for(n in 1:ncol(NormSuperList)){
  x <- NormSuperList[n]
  minMMidx <- which(x == min(x))
  
  # Finder 2nd lavest værdi hvis de to mindste ikke havde samme antal
  if( length(minMMidx < 2){
    min2MMidx <- which(x == min(x[x!=min(x)]))
    minMMidx[1:2] <- c(minMMidx,min2MMidx)
  }
  # Henter ref navne og sætter tabel op til tælling
  # Finder antal ref som har samme MM og ligger alle til tracker
  for(j in 1:length(minMMidx)){
    Ref_occur_tracker[minMMidx[j],1] <- Ref_occur_tracker[minMMidx[j],1] + 1
  }
}

min(NormSuperList[n[NormSuperList[n]!=min(NormSuperList[n])]])
x <- NormSuperList[n]







