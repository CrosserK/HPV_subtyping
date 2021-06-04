
# Script til at køre multiple fastq filer med HPV_subtyping.sh script
# Tager alle fastq filer der ligger i FastQF objekt og giver til HPV_subtyping.sh script en af gangen
# Kun "Define folders" sektion og kør. Vær sikker på at sigma_config_init.cfg ligger i Sigma_run mappe 

########## Define Folders ##########
SuperRunName=Run_len50_320
MainF=/home/pato/Skrivebord/HPV16_projekt
FastQF=$MainF/FASTQ
####################################


#Laver fastqnavne 
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
# Laver fastqnavne til RunName variabel i Sigma_run.sh
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

# Kører mega run (hardcoded antal filer ligenu)
END=$(wc -l $MainF/FASTQfiles_${SuperRunName}.txt)
START=1
for linenumber in {0..68} # $(eval echo "{$START..$END}")
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	HPV_subtyping.sh $FastqRunInput $FastqInput 
	echo $linenumber of 68 
done


######################################################################################################3


# Script til at køre multiple fastq filer med Specific_site_cov.sh script
# Tager alle fastq filer der ligger i FastQF objekt og giver til Specific_site_cov.sh script en af gangen
# Kun "Define folders" sektion og kør. Vær sikker på at sigma_config_init.cfg ligger i Sigma_run mappe 

########## Define Folders ##########
# FØR DENNE KØRES SKAL alle _filt fast filer fjernes
SuperRunName=Run_len50_320
MainF=/home/pato/Skrivebord/HPV16_projekt
FastQF=$MainF/FASTQ
####################################


#Laver fastqnavne 
#find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
# Laver fastqnavne til RunName variabel i Sigma_run.sh
#find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

# Kører mega run (hardcoded antal filer ligenu)
END=$(wc -l $MainF/FASTQfiles_${SuperRunName}.txt)
START=1
for linenumber in {1..68} # $(eval echo "{$START..$END}")
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	Specific_site_cov.sh $FastqRunInput $FastqInput
	echo $linenumber of 68
done