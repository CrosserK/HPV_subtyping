
# Call SNPs and create vcf file


##### EDIT #####
RunName=Pt_44.IonXpress_090_run
Reference=K02718.1
BamName=$Reference
#### EDIT DONE ####

MainF=/home/pato/Skrivebord/HPV16_projekt
SigmaF=$MainF/Sigma_run

Ref_FASTA=$MainF/References/$RunName/$Reference/${Reference}.fasta
BamFile=$SigmaF/$RunName/sigma_alignments_output/$BamName/Tester/${BamName}.align.bam
VarFile=$SigmaF/$RunName/sigma_alignments_output/$BamName/Tester/${BamName}.align_gatk.vcf

#freebayes -p 1 -f $Ref_FASTA $BamFile > $VarFile # -C for at bestemme hvor mange reads skal support en variant. -p for ploidy. Freebayes inkluderer leftalignment

#bcftools mpileup -f $Ref_FASTA $BamFile | bcftools call -v -Ov -c -o $VarFile

# Create sequence dictionary for gatk haplotypecaller
java -jar picard.jar CreateSequenceDictionary \
      R=$Ref_FASTA \
      O=${Ref_FASTA%.fasta}.dict

# Hvis 
java -jar picard.jar AddOrReplaceReadGroups \
I=$BamFile \
O=${BamFile%.bam}.fix.bam \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=ILY

mv $BamFile ${BamFile%.bam}.untouched.bam
mv ${BamFile%.bam}.fix.bam $BamFile

# Index bam
samtools index $BamFile

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $Ref_FASTA \
   -I $BamFile \
   -O $VarFile 
