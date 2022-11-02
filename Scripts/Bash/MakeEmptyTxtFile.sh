
# Laver tomme tekstfiler med samme navn som FASTQ filer i mappen FASTQ

MainF=$(dirname "$(readlink -f "$0")")
MainF=$MainF/../..

for f in $MainF/FASTQ/*.fastq; do
	Name=$(basename $f) 
	touch $MainF/References/InputRefs/${Name%fastq}txt
done