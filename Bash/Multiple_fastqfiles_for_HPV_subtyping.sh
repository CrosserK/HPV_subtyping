
# Script til at køre multiple fastq filer med HPV_subtyping.sh script
# Tager alle fastq filer der ligger i FastQF objekt og giver til HPV_subtyping.sh script en af gangen
# Kun "Define folders" sektion og kør. Vær sikker på at sigma_config_init.cfg ligger i Sigma_run mappe 

########## Define Folders ##########
SuperRunName=Exome_50_320_ampliconcalls_PaVE_revised
MainF=/home/pato/Skrivebord/HPV16_projekt
MinLen=50 # Min længde reads i fastqfil skal være
MaxLen=320 # max længde reads i fastqfil skal være
####################################

RunAnnoRSCript=true
RunSiteCovSCript=true
RunNoCallsSCript=true

FastQF=$MainF/FASTQ

#Laver fastqnavne 
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' > $MainF/FASTQfiles_${SuperRunName}.txt
# Laver fastqnavne til RunName variabel i Sigma_run.sh
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' > $MainF/FASTQfiles_${SuperRunName}_runnames.txt

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

############## ALIGNMENT AND VARIANT CALLING ##########################
START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	HPV_subtyping.sh $FastqRunInput $FastqInput $SuperRunName $MinLen $MaxLen
	echo $linenumber of $END
done
#######################################################################


# Følgende scripts læser hvilke FASTQfiler der skal behandles fra FASTQfiles_[navn].txt og FASTQfiles_[navn]_run.txt filerne

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder

Rscriptfolder=$MainF/Scripts/R
MultiFQFile=$(echo FASTQfiles_$SuperRunName)

# Finder navne for hver ref i References overmappe (eksl. undermapper)
find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_${SuperRunName}.txt
RefListFile=$RefdF/RefSubtyper_${SuperRunName}.txt
RefList=$(< $RefListFile)
SuperRunFolder=$SuperRunName


if [ $RunAnnoRSCript = true ]; then

	for refType in $RefList; do
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash.R $MainF/Annotation_results/ $refType $SuperRunName $MultiFQFile 
	done
fi

# Dette script skal bruge liste over alle referencer der køres

if [ $RunSiteCovSCript = true ]; then

	# Kører mega run
	START=1
	END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
	for (( linenumber=$START; linenumber<=$END; linenumber++ ));
	do
	
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	Specific_site_cov.sh $FastqRunInput $FastqInput $SuperRunFolder $RefListFile
	echo $linenumber of $END
	
	done

fi


# Kører R No_calls script 

if [ $RunNoCallsSCript = true ]; then

	for refType in $RefList; do
	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF/Annotation_results/ $refType $SuperRunName $MultiFQFile 
	# [SaveDir] [ReferenceName] [SuperRunName] [ListOfFastqFiles]
	done

fi


# Rscript $Rscriptfolder/Annotation_report.rmd $MainF/Annotation_results/ $SuperRunName $refType










############################################ Specific_site_cov.sh ##########################################################3


# Script til at køre multiple fastq filer med Specific_site_cov.sh script
# Tager alle fastq filer der ligger i FastQF objekt og giver til Specific_site_cov.sh script en af gangen
# Kun "Define folders" sektion og kør. Vær sikker på at sigma_config_init.cfg ligger i Sigma_run mappe 

########## Define Folders ##########
SuperRunName=Exome_50_120
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
