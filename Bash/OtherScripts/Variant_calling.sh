
# Call SNPs and create vcf file 


##### EDIT #####
RunName=Exome_50_320
Reference=K02718.1
FASTQName=Pt_153.IonXpress_003
BamName=K02718.1.sort.dup.readGroupFix.clipped.sort.sizefix.MD
#### EDIT DONE ####

MainF=/home/pato/Skrivebord/HPV16_projekt

Ref_FASTA=$MainF/References/${Reference}.fasta
BamFile=$MainF/Results/$RunName/${FASTQName}_run/$Reference/${BamName}.bam
VarFile=$MainF/Results/$RunName/${FASTQName}_run/$Reference/${BamName}_ploid1.vcf

#freebayes -p 1 -f $Ref_FASTA $BamFile > $VarFile # -C for at bestemme hvor mange reads skal support en variant. -p for ploidy. Freebayes inkluderer leftalignment

bcftools mpileup -f $Ref_FASTA $BamFile | bcftools call -v -Ov -c -o $VarFile

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
--sample-ploidy 1 \
   -R $Ref_FASTA \
   -I $BamFile \
   -O $VarFile \
   -bamout ${BamFile%.bam}_realigned.bam
