
# Script til at subtypere fastq filer i FASTQ mappe
# Kun ændre i "Indstillinger" sektion og kør enten gennem terminal eller copy paste ind i terminal 

########## Indstillinger ##########
# Put alle fastqfiler der skal testes ind i FASTQ mappe og sæt følgende:
Name=Karoline_92fastq  #Exome_50_320_ampliconcalls_PaVE_revised
MainF=/home/pato/Skrivebord/HPV16_projekt
VirStraindb=HPV16_16_virstrain_revised
QualTrim=20 # Qualtrim til cutadapt
MinLen=50 # Cutadapt fjerner alle reads under denne grænse
MaxLen=320 # Cutadapt fjerner alle reads over denne grænse
AmpliconRef=K02718.1 # Den reference som ampliconpanel er lavet ud fra. Muliggør en finere sortering af reads sammen med nedenstående bedfil som skal ufgøre amplicon position
BedFileNameX=IAD209923_226_Designed_compl # Bruges til at cutte primer sekvenser væk
####################################



########## Scripts to use ##########
indexReferences=false
cutOutsideAmplicons=false # Risk of loosing data if HPV is other type than the one used in ampliconpanel
qualityFilt=true # Skal ligenu være true, ellers skal fastqfiler gives endelsen "_filt.fastq" for at senere moduler kan finde dem
RunAnnoRSCript=true
RunSiteCovScript=true
RunNoCallsScript=true
####################################



Date=$(date +"%H%M_%d%m%Y")
SuperRunName=$(echo ${Name}_${Date})

Rscriptfolder=$MainF/Scripts/R
MultiFQFile=$(echo FASTQfiles_$SuperRunName)
SuperRunFolder=$SuperRunName

# Gå til MainF og opret mappe med samme navn som RunName, eks: 
mkdir -p $MainF/{Aligned/{Samfiles,BamFiles/unsorted},References/IndexedRef,Results/$SuperRunFolder/$RunName,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}} # -p option enables returning no error if folders exist. They will not be overwritten either way. 

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




############################# SUBTYPING ###############################
# Finder ligenu top 3 most possible strains
conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	VirStrain_subtyping.sh $FastqRunInput $FastqInput $SuperRunName $MainF $VirStraindb
	echo $linenumber of $END
done

conda deactivate
#######################################################################



####### REFERENCE INDEXERING TIL ALIGNMENT OG VARIANTCALLING ##########
# Putter hver reference i en undermappe for sig, så de er klar til kørsel
# genererer også dict fil og samtools faidx fil til HaplotypeCaller
# Genererer index til bwa mem
if [ $indexReferences = true ]; then

for f in $RefF/*.fasta; do
	RefName=${f##*/}
	mkdir -p $RefF/IndexedRef/${RefName%.fasta}
	cp $f $RefF/IndexedRef/${RefName%.fasta}/$RefName
	Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName

	# Create sequence dictionary for gatk haplotypecaller
	java -jar picard.jar CreateSequenceDictionary \
	      R=$Ref_FASTA \
	      O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt

	# Create index for bwa mem
	bwa index $Ref_FASTA

	# Index with samtools faidx
	samtools faidx $Ref_FASTA
done

fi
########################################################################



############## REVUDERING AF SUBTYPECALL, ALIGNMENT AND VARIANT CALLING ##########################
# Bruger _filt.fastq filer fra filtrering og SubtypeCall.txt fra subtypering. Selv hvis der ikke ønskes subtypering (så som hvis det er genotyperingspanel)
# der køres på, skal subtypering køres igennem for at lave disse filer. Derefter kan der manuelt vælges bestemte referencer inde i AlignAndVariantCall.sh scriptet.
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
	AlignAndVariantCall.sh $FastqRunInput $FastqInput $SuperRunName $MainF $VirRunOut_run $BedFileNameX $AmpliconRef
	echo fil $FastqInput færdig... $linenumber of $END
done
#######################################################################




# Følgende scripts læser hvilke FASTQfiler der skal behandles fra FASTQfiles_[navn].txt og FASTQfiles_[navn]_run.txt filerne

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
rm $MainF/FASTQ/*_filt.fastq;
#done



















######################### SIFT Amino acid substitution protein effect #########################
#OBS: HARDCODED TIL HPV16, søg HPV16-, hvis det skal ændres

## gff3 til protein fasta fil
#MainF=/home/pato/Skrivebord/HPV16_projekt
#
#File=$MainF/References/GFFfiles/K02718.1_revised.gff3
##Translationer:
#grep "CDS" $File | awk '{print $9}' | sed 's/.*translation=//' | sed 's/;.*//' > translations.txt
## tager alle linjer med CDS | printer 9. kolonne | tager alt efter translation= | tager alt før ";" 
## Protein navne:
#grep "CDS" $File | awk '{print $9}' | sed 's/.*protein_id=//' | sed 's/;.*//' > GeneNames.txt
#
## Samler til 1 fil
#
#
#rm -f $MainF/SIFT4G/Gene_fa.fasta
#touch $MainF/SIFT4G/Gene_fa.fasta
#START=1
#END=$(awk 'END{print NR}' translations.txt)
#for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
#	Gene=$(cat GeneNames.txt | awk -v var="$linenumber" 'NR==var {print}') 
#	Translation=$(cat translations.txt | awk -v var="$linenumber" 'NR==var {print}')
#	Header=$(echo \>${Gene})
#	echo -e "$Header""\n""$Translation" >> $MainF/SIFT4G/Gene_fa.fasta
#done
#
#rm -f translations.txt
#rm -f GeneNames.txt
#
#
## Oversætter fundne substitutioner til en protein fasta fil
##Fil med substitutioner
#FQFile=$MainF/Annotation_results/FASTQfiles_Karoline_run_test2_AnnotationFrequency_K02718.1_revised.txt
#sed 's/.*p.//' $FQFile | sed 's/ .*//' | awk 'NR!=1 {print}' > substitutions.txt
## Finder alt efter p. | finder alt før [space] | ignorerer første linje, da den er kolonne navn 
#cut -d ' ' -f3 $FQFile | awk 'NR!=1 {print}' > Genes.txt
#
#
## Generérer protein fasta fil for hver subs
#START=1
#END=$(awk 'END{print NR}' Genes.txt)
#
#for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
#	GeneToMod=$(awk -v var="$linenumber" 'NR==var {print}' Genes.txt)
#	SubPos=$(awk -v var="$linenumber" 'NR==var {print}' substitutions.txt | grep -o -E '[0-9]+')
#	# Finder linje | Finder position
#	SubAA=$(awk -v var="$linenumber" 'NR==var {print}' substitutions.txt | sed "s/.*$SubPos//")
#
#	# Tager header
#	grep $GeneToMod -A 1 $MainF/SIFT4G/Gene_fa.fasta | awk 'NR==1 {print}' > head.txt
#	# Finder gen transl og -A giver også næste N linjer | Fjerner første linje, da den er kolonne navn
#	grep $GeneToMod -A 1 $MainF/SIFT4G/Gene_fa.fasta | awk 'NR==2 {print}' | sed "s/./$SubAA/$SubPos" > ModdedGene.txt
#	cat head.txt ModdedGene.txt > $MainF/SIFT4G/Query_inst_${linenumber}.fasta
#	echo $(awk -v var="$linenumber" 'NR==var {print}' substitutions.txt ) >> $MainF/SIFT4G/SUBS/HPV16-${GeneToMod}.subst
#	rm -f head.txt
#	rm -f ModdedGene.txt
#done
#rm -f genes.txt
#rm -f substitutions.txt
#
#
## SIFT4G:
#for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
#linenumber=2
#/home/pato/sift4g/bin/sift4g -q $MainF/SIFT4G/Test_query.fasta --subst $MainF/SIFT4G/SUBS -d $MainF/SIFT4G/Gene_fa.fasta --out $MainF/SIFT4G/
## --outfmt light
#done


