
# Script til at køre multiple fastq filer med HPV_subtyping.sh script
# Tager alle fastq filer der ligger i FastQF objekt og giver til HPV_subtyping.sh script en af gangen
# Kun "Define folders" sektion og kør. Vær sikker på at sigma_config_init.cfg ligger i Sigma_run mappe 

########## Define Folders ##########
SuperRunName=Exome_50_120
MainF=/home/pato/Skrivebord/HPV16_projekt
MinLen=50 # Min længde reads i fastqfil skal være
MaxLen=120 # max længde reads i fastqfil skal være
####################################

FastQF=$MainF/FASTQ

#Laver fastqnavne 
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
# Laver fastqnavne til RunName variabel i Sigma_run.sh
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

# Kører mega run
START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	HPV_subtyping.sh $FastqRunInput $FastqInput $SuperRunName $MinLen $MaxLen
	echo $linenumber of $END
done






################################################ Specific_site_cov ######################################################3


# Script til at køre multiple fastq filer med Specific_site_cov.sh script
# Tager alle fastq filer der ligger i FastQF objekt og giver til Specific_site_cov.sh script en af gangen
# Kun "Define folders" sektion og kør. Vær sikker på at sigma_config_init.cfg ligger i Sigma_run mappe 

########## Define Folders ##########
# FØR DENNE KØRES SKAL alle _filt fast filer fjernes
SuperRunName=Exome_50_120
MainF=/home/pato/Skrivebord/HPV16_projekt
FastQF=$MainF/FASTQ
####################################


#Laver fastqnavne 
#find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
# Laver fastqnavne til RunName variabel i Sigma_run.sh
#find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	Specific_site_cov.sh $FastqRunInput $FastqInput $SuperRunName
	echo $linenumber of $END
done









##################### VIRSTRAIN ##############################
conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env

SuperRunName=qualtrim_50_320
MainF=/home/pato/Skrivebord/HPV16_projekt
FastQF=$MainF/FASTQ
VirStraindb=HPV16_16_virstrain # Skal være lavet med VirStrain_build.py
MinLen=50
MaxLen=120
####################################

# AKTIVER HVIS NYT RUN UDEN HPV_SUBTYPING.ŚH:
#Laver fastqnavne 
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
#Laver fastqnavne til RunName variabel i Sigma_run.sh
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

# Kører mega run (hardcoded antal filer ligenu)
START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	VirStrain_subtyping.sh $FastqRunInput $FastqInput $VirStraindb $MinLen $MaxLen
	echo $linenumber of $END
done

conda deactivate
