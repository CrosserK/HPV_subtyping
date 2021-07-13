

# ClIP TEST


MainF=/home/pato/Skrivebord/HPV16_projekt
Ref_FASTA=$MainF/References/K02718.1.fasta
BamFile=$MainF/Results/intersect_test/udenintersect.bam
BedFile=$MainF/References/BedFiles/IAD209923_226_Designed_poolx.bed


# Laver genome fil til bedtools complement
# awk -v OFS='\t' {'print $1,$2'} ${Ref_FASTA}.fai > ${Ref_FASTA%.fasta}_genome.txt
# bedtools complement -i $BedFile2 -g ${Ref_FASTA%.fasta}_genome.txt > ${BedFile2%.bed}_compl.bed

# Bruger amplicon regioner til at vælge reads
for i in 1 2; do

	# Choose file
	i=1
	tmpBedFile=${BedFile%x.bed}${i}.bed

	bedtools intersect -a $BamFile -b $tmpBedFile > ${BamFile%udenintersect.bam}intersect${i}.bam
	BamInt=${BamFile%udenintersect.bam}intersect${i}.bam
	#samtools ampliconclip --hard-clip -b ${tmpBedFile%.bed}_compl.bed $BamInt > ${BamInt%.bam}_clip.bam

	clipped_out=${BamInt}
	# TEST
	samtools sort $clipped_out > ${clipped_out%.bam}.sort.bam	

	samtools sort -n $clipped_out > ${clipped_out%.bam}.nsort.bam
	samtools fixmate -O bam ${clipped_out%.bam}.nsort.bam ${clipped_out%.bam}.nsort.sizefix.bam --verbosity 2
	samtools calmd ${clipped_out%.bam}.nsort.sizefix.bam $Ref_FASTA --output-fmt BAM > ${clipped_out%.bam}.nsort.sizefix.MD.bam
	samtools sort ${clipped_out%.bam}.nsort.sizefix.MD.bam > ${clipped_out%.bam}.sort.sizefix.MD.bam
	samtools index ${clipped_out%.bam}.sort.sizefix.MD.bam

	mv ${clipped_out%.bam}.sort.sizefix.MD.bam ${clipped_out%.bam}.sort.sizefix.MD_pool${i}.bam

done

# Samler splittede bamfiler
samtools merge ${clipped_out%.bam}.sort.sizefix.MD_pool1.bam ${clipped_out%.bam}.sort.sizefix.MD_pool2.bam









#### GAMMELT ¤#####








bedfile=/home/pato/Skrivebord/HPV16_projekt/References/IAD209923_226_Designed.bed
MainF=/home/pato/Skrivebord/HPV16_projekt/Results/Pt_33_RNA.IonXpress_087_run/K02718.1/
bam=$MainF/K02718.1.sort.dup.readGroupFix.bam
Ref_FASTA=/home/pato/Skrivebord/HPV16_projekt/References/K02718.1.fasta
clipped_out=$MainF/clipped_soft.bam






# Genome file til bedtools complement

#awk -v OFS='\t' {'print $1,$2'} ${Ref_FASTA}.fai > ${Ref_FASTA%.fasta}_genome.txt

#bedtools complement -i $bedfile -g ${Ref_FASTA%.fasta}_genome.txt > ${bedfile%.bed}_compl.bed

samtools ampliconclip --soft-clip --both-ends -b ${bedfile%.bed}_compl.bed $bam > $clipped_out
samtools sort -n $clipped_out > ${clipped_out%.bam}.nsort.bam
samtools fixmate -O bam ${clipped_out%.bam}.nsort.bam ${clipped_out%.bam}.nsort.sizefix.bam --verbosity 2
samtools calmd ${clipped_out%.bam}.nsort.sizefix.bam $Ref_FASTA --output-fmt BAM > ${clipped_out%.bam}.nsort.sizefix.MD.bam
samtools sort ${clipped_out%.bam}.nsort.sizefix.MD.bam > ${clipped_out%.bam}.sort.sizefix.MD.bam
samtools index ${clipped_out%.bam}.sort.sizefix.MD.bam