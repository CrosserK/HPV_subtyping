




# ClIP TEST


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