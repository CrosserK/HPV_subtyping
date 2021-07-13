
# Cut primer sequences from bed amplicon pool regions

RunName=$1
FastQFile=$2 #input uden extenstion
SuperRunFolder=$3
MainF=$4
Amplicon_ref=$5
BedFileNameX=$6
 
#Define Folders and params
BedFileName=$( echo  ${BedFileNameX}) # Regions in bedfile will be soft clipped from bam
BedFilePool=$MainF/References/BedFiles/${BedfileNameX}x.bed

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder

# Fetching Bed file
BedFile=$RefF/$BedFileName

# Genererer index til bwa mem
cp $RefF/${Amplicon_ref}.fasta $RefF/$RunName/${Amplicon_ref}/${Amplicon_ref}.fasta

Ref_FASTA=$RefF/$RunName/${Amplicon_ref}/${Amplicon_ref}.fasta

# Create index for bwa mem
bwa index $Ref_FASTA

# Index with samtools faidx
samtools faidx $Ref_FASTA