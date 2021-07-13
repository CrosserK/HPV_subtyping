
# Bed chrom name correcter
# Load bed fil uden header
bedfile <-read.table("/home/pato/Skrivebord/HPV16_projekt/References/IAD209923_226_Designed.bed")
bedfile[,1] <- "K02718.1"
write.table(bedfile, "/home/pato/Skrivebord/HPV16_projekt/References/IAD209923_226_Designedfix.bed", quote = F, col.names = F, row.names = F, sep = "\t")
