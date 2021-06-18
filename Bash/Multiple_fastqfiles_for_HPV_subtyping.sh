
# Script til at køre multiple fastq filer i FASTQ mappe
# Tager alle fastq filer der ligger i FastQF objekt og giver til HPV_subtyping.sh script en af gangen
# Kun ændre i "Define folders" sektion og kør 

########## Define Folders ##########
SuperRunName=KaroTest_rev3  #Exome_50_320_ampliconcalls_PaVE_revised
MainF=/home/pato/Skrivebord/HPV16_projekt
VirStraindb=HPV16_16_virstrain_revised
QualTrim=20 # Qualtrim til cutadapt
MinLen=50 # Min længde reads i fastqfil skal være
MaxLen=320 # max længde reads i fastqfil skal være
AmpliconRef=K02718.1
BedFileNameX=Revised_IAD209923_226_Designed_compl # Den bed fil som amplicons er lavet ud fra i vådlab del. Bruges til at cutte primer sekvenser væk
####################################

cutOutsideAmplicons=false # Risk of loosing data if HPV is other type than the one used in ampliconpanel
RunAnnoRSCript=true
RunSiteCovSCript=true
RunNoCallsSCript=true

# Gå til MainF og opret mappe med samme navn som RunName, eks: 
mkdir -p $MainF/{Aligned/{Samfiles,BamFiles/unsorted},References,Results/$SuperRunFolder/$RunName,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}} # -p option enables returning no error if folders exist. They will not be overwritten either way. 

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder

FastQF=$MainF/FASTQ
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' \
> $MainF/FASTQfiles_${SuperRunName}.txt #Finder FASTQ filer i FASTQ mappe og laver tekstfil med alle navne 
find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' \
> $MainF/FASTQfiles_${SuperRunName}_runnames.txt # Laver fastqnavne til RunName variabel  (så den er anderledes en navnet)

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) # Henter liste og gemmer som variabel
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

########################## Primer sekvens cut #########################
if [ $cutOutsideAmplicons = true ]; then
	
	START=1
	END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 	
	PrimerSeqCut.sh
 	done
fi




############################# SUBTYPING ###############################
conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	VirStrain_subtyping.sh $FastqRunInput $FastqInput $SuperRunName $MainF $VirStraindb $QualTrim $MinLen $MaxLen
	echo $linenumber of $END
done

conda deactivate
#######################################################################


####### REFERENCE INDEXERING TIL ALIGNMENT OG VARIANTCALLING ##########
# Putter hver reference i en undermappe for sig, så de er klar til kørsel
# genererer også dict fil og samtools faidx fil til HaplotypeCaller
# Genererer index til bwa mem
for f in $RefF/*.fasta; do
	mkdir -p ${f%.fasta}
	RefName=${f##*/}
	cp $f IndexedRef/${f%.fasta}/$RefName
	Ref_FASTA=IndexedRef/${f%.fasta}/$RefName

	# Create sequence dictionary for gatk haplotypecaller
	java -jar picard.jar CreateSequenceDictionary \
	      R=$Ref_FASTA \
	      O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt

	# Create index for bwa mem
	bwa index $Ref_FASTA

	# Index with samtools faidx
	samtools faidx $Ref_FASTA
done
########################################################################


# Bruger _filt.fastq filer fra SUBTYPING
############## ALIGNMENT AND VARIANT CALLING ##########################
START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	VirSupOut=$MainF/VirStrain_run/$SuperRunName
	RunName=${FastqInput}_run
	VirRunOut_run=$VirSupOut/$RunName
	AlignAndVariantCall.sh $FastqRunInput $FastqInput $SuperRunName $MainF $VirRunOut_run $BedFileNameX $AmpliconRef
	echo fil $FastqInput færdig... $linenumber of $END
done
#######################################################################

# Ovenstående returnerer RevRefCalls variabel, som er den korrekte fundne subtyper

# Følgende scripts læser hvilke FASTQfiler der skal behandles fra FASTQfiles_[navn].txt og FASTQfiles_[navn]_run.txt filerne

Rscriptfolder=$MainF/Scripts/R
MultiFQFile=$(echo FASTQfiles_$SuperRunName)

# Finder navne for hver ref i References overmappe (eksl. undermapper)
#find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
#sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_${SuperRunName}.txt
#RefListFile=$RefdF/RefSubtyper_${SuperRunName}.txt
RefList=$RevRefCalls # Henter variabel fra alignement/variantcalling script    # $(< $RefListFile)
SuperRunFolder=$SuperRunName


if [ $RunAnnoRSCript = true ]; then
	START=1
	END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
	for (( linenumber=$START; linenumber<=$END; linenumber++ )); 
	do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	VirSupOut=$MainF/VirStrain_run/$SuperRunName
	RunName=${FastqInput}_run
	VirRunOut_run=$VirSupOut/$RunName
	RevRefCalls=$(< $VirRunOut_run/Revised_SubTypeCalls.txt) 
	for refType in $RevRefCalls; do
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash.R $MainF/Annotation_results/ $refType $SuperRunName $MultiFQFile
	done
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
















