
# Laver tomme tekstfiler med samme navn som FASTQ filer i mappen FASTQ

MainF=/home/pato/Skrivebord/HPV16_projekt

for f in $MainF/FASTQ/*.fastq; do
	Name=$(basename $f) 
	touch $MainF/References/InputRefs/${Name%fastq}txt
done