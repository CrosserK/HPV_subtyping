
# Cut primer sequences from bed amplicon pool regions

set -o errexit # Exit if something fails
set -u nounset # Exit if undeclared variables

# Reading input parameters
while [[ "$#" -gt 0 ]]; do
	case "$1" in
        -s|--superrunname) SuperRunName="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -m|--mainf) MainF="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -a|--ampl) Amplicon_ref="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -b|--bed) BedFileNameX="$2"; [ "$(echo "$2" | cut -c 1)" == "-" ] || [ "$2" == "" ] && help || shift 2;;
        -h|--help) help;;
        *) help;;
    esac
done

#Define Folders and params
BedFileName=$( echo  ${BedFileNameX}) # Regions in bedfile will be soft clipped from bam
BedFilePool=$MainF/References/BedFiles/${BedfileNameX}x.bed
RefF=$MainF/References

# Fetching Bed file
BedFile=$RefF/$BedFileName

# Genererer index til bwa mem
cp $RefF/${Amplicon_ref}.fasta $RefF/$SuperRunName/${Amplicon_ref}/${Amplicon_ref}.fasta

Ref_FASTA=$RefF/$SuperRunName/${Amplicon_ref}/${Amplicon_ref}.fasta

# Create index for bwa mem
bwa index $Ref_FASTA

# Index with samtools faidx
samtools faidx $Ref_FASTA