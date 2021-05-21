
#Target folders:
TargetFolder='/home/pato/Skrivebord/HPV16_projekt'
RefFolder=$TargetFolder/References

cd $RefFolder
cat HPV16_multi.fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,19,11) ".fasta")}
        print $0 > filename
}'


