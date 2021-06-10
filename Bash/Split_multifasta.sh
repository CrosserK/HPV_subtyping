
#Target folders:
TargetFolder='/home/pato/Skrivebord/HPV16_projekt'
RefFolder=$TargetFolder/References/Bruges_ikke

cd $RefFolder
cat 15_HPV_sublineages.fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2,11) ".fasta")}
        print $0 > filename
}'


