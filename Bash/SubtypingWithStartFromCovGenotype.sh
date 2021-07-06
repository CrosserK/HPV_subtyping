
# Genererer custom referencer efter genotype calls, hvis der er flere HPV typer i en prøve

###############
Name=Karoline_run_covGenTest 
MainF=/home/pato/Skrivebord/HPV16_projekt
QualTrim=20 # Qualtrim til cutadapt
MinLen=50 # Cutadapt fjerner alle reads under denne grænse
MaxLen=320 # Cutadapt fjerner alle reads over denne grænse
AmpliconRef=K02718.1 # Den reference som ampliconpanel er lavet ud fra. Muliggør en finere sortering af reads sammen med nedenstående bedfil som skal ufgøre amplicon position
BedFileNameX=IAD209923_226_Designed_compl # Bruges til at cutte primer sekvenser væk
################


########## Scripts to use ##########
indexReferences=false
cutOutsideAmplicons=false # Risk of loosing data if HPV is other type than the one used in ampliconpanel
qualityFilt=true # Skal ligenu være true, ellers skal fastqfiler gives endelsen "_filt.fastq" for at senere moduler kan finde dem
RunAnnoRSCript=true
RunSiteCovScript=true
RunNoCallsScript=true
####################################

# Date= # Uncomment til testfiler og comment næste linje
Date=$(date +"%H%M_%d%m%Y")
SuperRunName=$(echo ${Name}_${Date})

Rscriptfolder=$MainF/Scripts/R
MultiFQFile=$(echo FASTQfiles_$SuperRunName)

####################################

# Gå til MainF og opret mappe med samme navn som RunName, eks: 
mkdir -p $MainF/{Aligned/{Samfiles,BamFiles/unsorted},GenotypeCalls,References/{IndexedRef,Combined_refs},Results/$SuperRunFolder/$RunName,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}} # -p option enables returning no error if folders exist. They will not be overwritten either way. 

SeqF=$MainF/FASTQ; AnaF=$MainF/Analysis; QualF=$AnaF/Qual; DepthF=$AnaF/Depth; FlagF=$AnaF/Flagstats;
DupF=$AnaF/DuplicateMetrics; RefdF=$MainF/ReferenceDetails; RefF=$MainF/References; ErrorF=$AnaF/Errors; ResultsF=$MainF/Results/$SuperRunFolder; InputRefs=$SeqF/InputRefsTest

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

############################# FASTQ FILTRERING ###############################
if [ $qualityFilt = true ]; then

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do

	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	FastqInputAddr=$SeqF/${FastqInput}.fastq
	cutadapt -q $QualTrim -m $MinLen -M $MaxLen ${FastqInputAddr} -o ${FastqInputAddr%.fastq}_filt.fastq
	
	done

fi
# Dette navngiver filer med endelsen "_filt.fastq", hvilket er et krav for at resten af workflowet kan finde fastq filerne
##############################################################################


############## LAVER KOMBINEREDE REFERENCER #########################

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
for (( linenumber=$START; linenumber<=$END; linenumber++ )); do

	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	FastqInputAddr=$SeqF/${FastqInput}.fastq

	# Finder antal referencer i referencefil
	RefNumber=$(wc -l $InputRefs/${FastqInput}.txt | awk '{ print $1 }')
	
	# Hvis mere end 1 type i prøve:

	if [ $RefNumber -gt 1 ]; then

	# Læser referencer som array
	readarray -t RefsArr < $InputRefs/${FastqInput}.txt

	# Initerer array
	declare -a AllAddressList="" 

	rm -f $MainF/GenotypeCalls/${FastqInput}_SplitTo.txt
	for ref in ${RefsArr[@]}; do
	
	AllAddressList+=" "
	AllAddressList+=$( find $RefF -maxdepth 1  -name "${ref}*" -type f )
	echo "${ref}" >> $MainF/GenotypeCalls/${FastqInput}_SplitTo.txt

	done

	# Laver nyt fasta navn for kombineret reference
	newRefName=$(printf "_%s" "${RefsArr[@]}")
	newRefName=${newRefName:1} # Fjerner forkert foranstående "_" der opstår i ovenstående
	newRefNameFasta=${newRefName}.fasta

	# Send refname til fundne referencer for fil 
	echo $newRefName > $MainF/GenotypeCalls/${FastqInput}.txt

	# Har nu adresser til pågældende referencer. Samler:
	cat $AllAddressList > $RefF/Combined_refs/$newRefNameFasta


	# INDEXERER
	NewRef=$RefF/Combined_refs/$newRefNameFasta
	RefName=${NewRef##*/}
	mkdir -p $RefF/IndexedRef/${RefName%.fasta}
	cp $NewRef $RefF/IndexedRef/${RefName%.fasta}/$RefName
	Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName

	# Create sequence dictionary for gatk haplotypecaller
	java -jar picard.jar CreateSequenceDictionary \
	      R=$Ref_FASTA \
	      O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt

	# Create index for bwa mem
	bwa index $Ref_FASTA

	# Index with samtools faidx
	samtools faidx $Ref_FASTA



	else
	# Ellers hvis der kun er en type i fastq fil, gem da referencenavn i variabel
	Ref_FASTA=$(< $InputRefs/${FastqInput}.txt)
	echo "${ref}" > $MainF/GenotypeCalls/${FastqInput}_SplitTo.txt

	fi 

done


################## MULTIPLE TYPES AWARE ALIGNMENT & VARIANT CALLING ###########################
# Bruger _filt.fastq filer fra filtrering og referencekald fra KOMBINEREDE REFERENCER modul. 
# Laver 1 fil med alle referencer til No_call_script
find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_latest.txt

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	VirSupOut=$MainF/VirStrain_run/$SuperRunName
	RunName=${FastqInput}_run
	VirRunOut_run=$VirSupOut/$RunName
	AlignAndVariantCallFromCovGenotyping.sh $FastqRunInput $FastqInput $SuperRunName $MainF $VirRunOut_run $BedFileNameX $AmpliconRef
	echo fil $FastqInput færdig... $linenumber of $END
done
################################################################################################




# Split i hver funden type 

# FOR EACH UNIQUE CHRNAME SPLIT

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do



done

#######################################################################




# Variantcall







# Følgende scripts læser hvilke FASTQfiler der skal behandles fra FASTQfiles_[navn].txt og FASTQfiles_[navn]_run.txt filerne

# Rscript -e "library(GenomicFeatures)"


if [ $RunAnnoRSCript = true ]; then

	# Dette script tager selv fat i nødvendige Fastfiler fra MultiFQFile
	VirSupOut=$MainF/VirStrain_run/$SuperRunName
	RunName=${FastqInput}_run
	VirRunOut_run=$VirSupOut/$RunName
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_v2.R $MainF $MainF/Annotation_results $SuperRunName $MultiFQFile

fi

# Dette script skal bruge liste over alle referencer der køres. Fungerer kun hvis alle fastq filer har været alignet til samme referencer

if [ $RunSiteCovScript = true ]; then

	START=1
	END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)
	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do

	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	Specific_site_cov.sh $FastqRunInput $FastqInput $SuperRunFolder 
	
	done

fi


# Kører R No_calls script 

if [ $RunNoCallsScript = true ]; then

	#for refType in $RefList; do

	#if [ -d $ResultsF/$FastqRunInput/$refType ]; then
	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF $MainF/Annotation_results $SuperRunName $MultiFQFile
	# [SaveDir] [ReferenceName] [SuperRunName] [ListOfFastqFiles]
	
	#fi

	#done
	#done

fi

# Fjerner alle filtrerede fastqfiler der blev genereret, så de ikke bliver filtreret en gang til ved ny kørsel og der opstår dobbeltfiler
# rm $MainF/FASTQ/*_filt.fastq;
# done