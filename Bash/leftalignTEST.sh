while getopts u:a:f: flag
do
    case "${flag}" in
        u) username=${OPTARG};;
        a) age=${OPTARG};;
        f) fullname=${OPTARG};;
    esac
done
echo "Username: $username";
echo "Age: $age";
echo "Full Name: $fullname";

Ref_FASTA='/home/pato/Skrivebord/HPV16_projekt/References/runLeftAlTest2/K02718.1/K02718.1.fasta' 
insam='/home/pato/Skrivebord/HPV16_projekt/Sigma_run/runLeftAlTest2/sigma_alignments_output/K02718.1/K02718.1.align.sam' 

outbam='/home/pato/Skrivebord/HPV16_projekt/Sigma_run/runLeftAlTest2/sigma_alignments_output/K02718.1/K02718.1.align.sort.bam' 


# Sort bamfile
samtools sort $insam -o $outbam
	
gatk LeftAlignIndels \
  -R $Ref_FASTA \
  -I $outbam \
  -O ${outbam%bam}leftAl1.bam \
  --disable-tool-default-read-filters