

#Define Folders
SuperRunName=Megakørsel1 
MainF=/home/pato/Skrivebord/HPV16_projekt
##########


#Laver fastqnavne 
FastQF=$MainF/FASTQ
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
# Laver fastqnavne til RunName variabel i Sigma_run.sh
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

# Kører mega run (hardcoded antal filer ligenu)
END=$(wc -l $MainF/FASTQfiles_${SuperRunName}.txt)
START=1
for linenumber in {1..68} # $(eval echo "{$START..$END}")
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	Sigma_run.sh $FastqRunInput $FastqInput
done