
# INDEXERER OG LÆGGER I IndexedRefs MAPPE

FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre
RefF=/home/pato/Skrivebord/HPV16_projekt/References

for f in $FilePath/HPV16_16subtypes/*.fasta; do

	RefName=${f##*/}
	echo $RefName

	mkdir -p $RefF/IndexedRef/${RefName%.fasta}
	cp $f $RefF/IndexedRef/${RefName%.fasta}/$RefName
	Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName

	# Create sequence dictionary for gatk haplotypecaller
	java -jar picard.jar CreateSequenceDictionary \
	R=$Ref_FASTA \
	O=${Ref_FASTA%.fasta}.dict

	# Create index for bwa mem
	bwa index $Ref_FASTA

	# Index with samtools faidx
	samtools faidx $Ref_FASTA

done