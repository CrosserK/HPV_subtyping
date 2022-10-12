#!/bin/bash

# HPV_subtyping main script, that controls and calls all necessary functions

function parse_yaml {

	# Enables yaml file to be turned into objects

   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# Read all yaml objects into environment
eval $(parse_yaml Configurations.yaml)

set -u  # Exit if unset variable is called, must be after parse_yaml function.
set -e  # Exit if something exits with a non-zero status

DIR_SCRIPT_BASE=$(dirname "$(readlink -f "$0")")

# Setting SuperName for run
if [ $customName = true ]; then
	SuperRunName=$(grep "cName=" <<< $Options | sed 's/cName=//g') 
	echo "Using custom name and saved reference info..."
else
	# append date and time to name
	Date=$(date +"%H%M_%d%m%Y")
	SuperRunName=$(echo ${RunName}_${Date})
	FastQF=$MainF/FASTQ
	# Listing fastq names in a file
	find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' \
	> $MainF/FASTQfiles_${SuperRunName}.txt    # Finding FASTQ files in FASTQ folder and create text file with them all
	find $FastQF/ -maxdepth 1  -name '*.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' \
	> $MainF/FASTQfiles_${SuperRunName}_runnames.txt
fi

# Gå til MainF og opret mappe med samme navn som RunName: 
mkdir -p $MainF/{GenotypeCalls,References/{IndexedRef,Combined_refs},Results/$SuperRunName/,ReferenceDetails,FASTQ,Analysis/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}}  

# Defining folder locations
RunFiles=$MainF/FASTQfiles_${SuperRunName}.txt
RunFilesR=$MainF/FASTQfiles_${SuperRunName}_runnames.txt
SeqF=$MainF/FASTQ; 
AnaF=$MainF/Analysis; 
RefdF=$MainF/ReferenceDetails; 
RefF=$MainF/References; 
ErrorF=$AnaF/Errors; 
InputRefs=$RefF/InputRefs;
Rscriptfolder=$MainF/Scripts/R;
BashScriptF=$MainF/Scripts/Bash;
PythonScriptF=$MainF/Scripts/Python;
MultiFQFile=$(echo FASTQfiles_$SuperRunName);
ResultsF=$MainF/Results/$SuperRunName;

############# Henter fastq liste og gemmer som variabel ##############
FastqList=$(< $RunFiles) 
FastqRunList=$(< $RunFilesR)

################# Finder genotype fra covarage matrix #################
rm -f $RefF/*fasta.fai
if [ $genotypeFromCovMatrix = true ]; then
	python $PythonScriptF/autoDetectTypeFromCov.py -i $SeqF/covMatrix.csv -o $InputRefs -r $RefF
fi

########################## Primer sequence cut (BETA) #########################
if [ $cutOutsideAmplicons = true ]; then
	START=1
	END=$(awk 'END{print NR}' $RunFiles)
	# One by one, cut primer sequences
	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 	
	"${BashScriptF}"/PrimerSeqCut.sh -s "${SuperRunName}" -m "$MainF" -a "$AmpliconRef" -b "$BedFileNameX"  
 	done
fi

############################# FASTQ FILTRERING ###############################
if [ $qualityFilt = true ]; then
START=1
END=$(awk 'END{print NR}' $RunFiles)
# One by one, filter with cutadapt
for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
	FastqInputAddr=$SeqF/${FastqInput}.fastq
	# This gives filtered files the ending _filt.fastq
	cutadapt -q $QualTrim -m $MinLen -M $MaxLen "${FastqInputAddr}" -o "${FastqInputAddr%.fastq}"_filt.fastq
done
	# Update the list of files with filt names
	find $FastQF/ -maxdepth 1  -name '*_filt.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' \
	> $MainF/FASTQfiles_${SuperRunName}.txt    # Finding FASTQ files in FASTQ folder and create text file with them all
	find $FastQF/ -maxdepth 1  -name '*_filt.fastq' | sed 's/^.*FASTQ.//' | sed 's/.fastq//' | sed 's/$/_run/' \
	> $MainF/FASTQfiles_${SuperRunName}_runnames.txt
echo done filtering fastq files
fi

####### REFERENCE INDEXING FOR ALIGNMENT AND VARIANTCALLING ##########
# Put each reference in a subfolder and generates variuos dict files needed
if [ $indexReferences = true ]; then
for f in $RefF/*.fasta; do
	RefName=${f##*/}
	mkdir -p $RefF/IndexedRef/${RefName%.fasta}
	# Check if fasta is already indexed
	if [ -f $RefF/IndexedRef/${RefName%.fasta}/${RefName} ]; then
	echo $RefName already in IndexedRef folder
	else
	cp $f $RefF/IndexedRef/${RefName%.fasta}/$RefName
	Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName
	samtools faidx $Ref_FASTA
	# Create sequence dictionary for gatk haplotypecaller
	java -jar ~/picard.jar CreateSequenceDictionary \
	R=$Ref_FASTA \
	O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt
	bwa index $Ref_FASTA
	fi
done
fi

# Setting input references for all fastq files if they are set as the same
if [ $customRefForAll = true ]; then
	START=1
	END=$(awk 'END{print NR}' $RunFiles)
	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
	rm -f $InputRefs/${FastqInput%_filt}.txt
	touch $InputRefs/${FastqInput%_filt}.txt
	for newref in ${cRef[@]}; do
	echo $newref >> $InputRefs/${FastqInput%_filt}.txt
	done
	done
fi


################ VIRSTRAIN SUBTYPING  #################
# (genotyped with covMatrix or manually entered in text files in References/InputRefs)
if [ $VirStrainSubTyping = true ]; then
mkdir -p $MainF/GenotypeCalls/$SuperRunName/
START=1
END=$(awk 'END{print NR}' $RunFiles)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR)
	FastqInputAddr=$SeqF/${FastqInput}.fastq 
	#  Getting correct VirStrain subtypedatabase for genotype
	mkdir -p $MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}
	outFolder=$MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}
	MainCallFile=$InputRefs/${FastqInput%_filt}.txt
	#  For each genotype in inputref file, find subtype
	genoStart=1
	genoEnd=$(awk 'END{print NR}' "$InputRefs"/"${FastqInput%_filt}".txt)
	for (( genonum=$genoStart; genonum<=$genoEnd; genonum++ )); do
	VirStrainSub=$(sed -n "${genonum}"p "$MainCallFile" | grep -o "^HPV[0-9]*" ) # -m to only search defined line and -o to only output match
	MainCallFullName=$(sed -n "${genonum}"p "$MainCallFile" | grep "^HPV[0-9]*")
	echo Finding $VirStrainSub subtype for ${FastqInput}...
	BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
	VirStrainSubDB=$MainF/References/$VirStrainSubFolders/$BaseDBName
	"${BashScriptF}"/VirStrain_genotypeToSubtype_covCalled.sh -r $FastqRunInput -f $FastqInputAddr -s $SuperRunName -m $MainF -n $genonum -v $VirStrainSubDB -c $MainCallFullName
	done
done
fi

# ############################# VIRSTRAIN GENOTYPING & SUBTYPING ###############################

if [ $VirStrainGenoAndSubTyping = true ]; then
# Finding top 3 most possible strains, using top 1 for alignment
VSdb=$MainF/References/$VirStrainMaindb
START=1
END=$(awk 'END{print NR}' $RunFiles)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
	"${BashScriptF}"/VirStrain_genotyping.sh $FastqRunInput $FastqInput $SuperRunName $MainF $VSdb
	echo $linenumber of $END
done;

# SUBTYPING
START=1
END=$(awk 'END{print NR}' $RunFiles)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
	#  Getting correct VirStrain subtypedatabase for genotype
	SuperRunOut=$MainF/VirStrain_run/$SuperRunName/GenotypeCalls
	Resultsout=$SuperRunOut/$FastqRunInput
	MainCallFile=$Resultsout/SubtypeCall.txt
	VirStrainSub=$(grep -m1 -o "^HPV[0-9]*" $MainCallFile) # -m to only search defined line and -o to only output match
	MainCallFullName=$(grep -m1 "^HPV[0-9]*" $MainCallFile)
	BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
	VirStrainSubDB=$MainF/References/$VirStrainSubFolders/$BaseDBName
	"${BashScriptF}"/VirStrain_genotypeToSubtype.sh -r $FastqRunInput -f $FastqInputAddr -s $SuperRunName -m $MainF -v $VirStrainSubDB -c $MainCallFullName
	echo $linenumber of $END
done
fi


############## Combining references and/or saving in genotypcalls folder #################
mkdir -p $MainF/GenotypeCalls/$SuperRunName
if [ $CombineRefs = true ] && [ $VirStrainGenoAndSubTyping = false ]; then
START=1
END=$(awk 'END{print NR}' $RunFiles)
for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
FastqInput=$(awk "NR==$linenumber" $RunFiles)
FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
FastqInputAddr=$SeqF/${FastqInput}.fastq
# Finding number of references in reference file
RefNumber=$(wc -l $MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}.txt | awk '{ print $1 }')
# If more than one type in sample, combine reference files (fasta) and save info in GenotypeCalls folder
if [ $RefNumber -gt 1 ]; then
	# Read references as array
	readarray -t RefsArr < $MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}.txt
	# Initiate array
	declare -a AllAddressList="" 
	for ref in ${RefsArr[@]}; do
	# Find ref path and append to list
	AllAddressList+=" "
	AllAddressList+=$( find $RefF -maxdepth 1  -name "${ref}*" -type f )
	#echo "${ref}" >> $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}_SplitTo.txt
	done
	# Make new fasta name for combined ref
	newRefName=$(printf "_%s" "${RefsArr[@]}")
	newRefName=${newRefName:1} # Removes wrong prefix of "_" from above method
	newRefNameFasta=${newRefName}.fasta
	# Collect the fastas to one file
	cat $AllAddressList > $RefF/Combined_refs/$newRefNameFasta
	# Index and save in IndexedRefs folder
	NewRef=$RefF/Combined_refs/$newRefNameFasta
	RefName=${NewRef##*/}
	mkdir -p $RefF/IndexedRef/${RefName%.fasta}
	cp $NewRef $RefF/IndexedRef/${RefName%.fasta}/$RefName
	Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName
	# Create sequence dictionary for gatk haplotypecaller
	java -jar ~/picard.jar CreateSequenceDictionary \
	R=$Ref_FASTA \
	O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${RefName}.txt
	# Create index for bwa mem
	bwa index $Ref_FASTA
	# Index with samtools faidx
	samtools faidx $Ref_FASTA
else
	# Else if only one detected hpvtype in fastq, save reference from txt file in variable
	Ref_FASTA=$(< $MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}.txt)
fi 
done

# elif [ $CombineRefs = false ] && [ $VirStrainGenoAndSubTyping = false ] && [ $customName = false ]; then

# 	# Ellers hvis der ikke skal kombineres referencer, gem da referencenavne fra inputrefs i genotypecalls
# 	START=1
# 	END=$(awk 'END{print NR}' $RunFiles)
# 	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
# 	FastqInput=$(awk "NR==$linenumber" $RunFiles)
# 	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 

# 	rm -f $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}.txt
# 	touch $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}.txt

# 	rm -f $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}_SplitTo.txt
# 	touch $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}_SplitTo.txt

# 	for newref in ${cRef[@]}; do
# 	echo $newref >> $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}.txt
# 	echo $newref >> $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}_SplitTo.txt
# 	done
# 	done
fi
###############################################################################################
# Output: $MainF/GenotypeCalls/${SuperRunName}/${FastqInput}.txt and _SplitTo.txt version

################## MULTIPLE TYPES AWARE ALIGNMENT & VARIANT CALLING ###########################
# Input: References from $MainF/GenotypeCalls/$SuperRunName/${FastQFile}.txt and fastq files 

if [ $AlignAndVarCall = true ]; then 
START=1
END=$(awk 'END{print NR}' $RunFiles)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
	FastqInputAddr=$SeqF/${FastqInput}.fastq 
	# TEST
	# GET EACH TYPE FOUND FOR EACH FASTQ !!
	GenoTypeCallIn=$MainF/GenotypeCalls/$SuperRunName
	
	"${BashScriptF}"/AlignAndVariantCall.sh -r $FastqRunInput -f $FastqInputAddr -s $SuperRunName -m $MainF -g $GenoTypeCallIn -b $BedFileNameX -a $AmpliconRef
	echo File $FastqInput done with alignment and variant calling... $linenumber of $END fastq files
done
fi
# Fjerner alle filtrerede fastqfiler der blev genereret, så de ikke bliver filtreret en gang til ved ny kørsel og der opstår dobbeltfiler
rm $MainF/FASTQ/*_filt.fastq $RunFiles
################################################################################################

# Moving files for convenience, if subtyping has been done
if [ $VirStrainGenoAndSubTyping = true ]; then
	FILE1=$MainF/VirStrain_run/$SuperRunName/VirStrain_summary.txt
	FILE2=$MainF/VirStrain_run/$SuperRunName/GenotypeCalls/VirStrain_summary.txt
	cp $FILE1 $MainF/Results/$SuperRunName/VirStrainSubtypeCalls.txt
	cp $FILE2 $MainF/Results/$SuperRunName/VirStrainGenotypeCalls.txt
fi

# # FOR EACH UNIQUE CHRNAME: SPLIT & VARIANTCALL AGAIN
# if [ $CombineRefs = true ]; then
# START=1
# END=$(awk 'END{print NR}' $RunFiles)

# for (( linenumber=$START; linenumber<=$END; linenumber++ ))
# do
# 	FastqInput=$(awk "NR==$linenumber" $RunFiles)
# 	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 

# 	# Now combine maintypes to split reads into different maintypes
# 	GenoTypeCallIn=$MainF/GenotypeCalls/$SuperRunName
# 	GenoTypeCallForFastq=$(< $MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}.txt)	
	
# 	BamName=${FastqInput}_${GenoTypeCallForFastq}.sort.bam
# 	bamtools split -in $ResultsF/$FastqRunInput/$GenoTypeCallForFastq/ResultFiles/$BamName -reference
	
# 	#  Samler hver splittet fil i en mappe og laver variantcalls
# 	StartRefs=1
# 	EndRefs=$(awk 'END{print NR}' $MainF/GenotypeCalls/$SuperRunName/${FastqInput}_SplitTo.txt)
# 	if [ $EndRefs -gt 1 ]; then 
# 	for (( linen=$StartRefs; linen<=$EndRefs; linen++ )); do
# 	SplitRef=$(awk "NR==$linen" $MainF/GenotypeCalls/$SuperRunName/${FastqInput}_SplitTo.txt)
# 	mkdir -p $ResultsF/$FastqRunInput/$SplitRef
# 	SplitBam=$ResultsF/$FastqRunInput/$GenoTypeCallForFastq/ResultFiles/${FastqInput}_${GenoTypeCallForFastq}.sort.REF_${SplitRef}.bam
# 	mv $SplitBam $ResultsF/$FastqRunInput/$SplitRef/${FastqInput}_${SplitRef}.bam

# 	#  Variant filtrering, indexering og sortering
# 	"${BashScriptF}"/VariantCallFromCovGenotypingSplitRefs.sh $FastqRunInput $FastqInput $SuperRunName $MainF $GenoTypeCallIn $SplitRef $BedFileNameX $AmpliconRef 
# 	done
# 	fi
# 	echo fil $FastqInput færdig med split calls... $linenumber af $END fastqfiler
# done
# fi
#######################################################################


# Fixer navn på vcf filer som har fået fjernet alle filtrerede varianter helt, hvis der ikke har været kombinerede referencer med splittede bamfiler
START=1
END=$(awk 'END{print NR}' $RunFiles)
for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 

	FILE1=$MainF/GenotypeCalls/$SuperRunName/${FastqInput}_SplitTo.txt
	FILE2=$MainF/GenotypeCalls/$SuperRunName/${FastqInput%_filt}.txt

	# Checker om en fastqfil havde kombinerede referencer eller ej
	if cmp --silent -- "$FILE1" "$FILE2"; then
	StartRefs=1
	EndRefs=$(awk 'END{print NR}' $MainF/GenotypeCalls/$SuperRunName/${FastqInput}_SplitTo.txt)
	for (( linen=$StartRefs; linen<=$EndRefs; linen++ )); do
	#linen=1
	SplitRef=$(awk "NR==$linen" $MainF/GenotypeCalls/$SuperRunName/${FastqInput}_SplitTo.txt)
	# renamer så rscript kan finde vcf fil
	#SplitRef=$(awk "NR==1" $MainF/GenotypeCalls/$SuperRunName/${FastqInput}_SplitTo.txt)
	# Her indsættes hvis der skal bruge dup markerede ".dup."
	cp $ResultsF/$FastqRunInput/${SplitRef}/${SplitRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf $ResultsF/$FastqRunInput/${SplitRef}/${FastqInput}_${SplitRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf
	done
	fi

done


#  Følgende scripts læser hvilke FASTQfiler der skal behandles fra FASTQfiles_[navn].txt og FASTQfiles_[navn]_run.txt filerne
if [ $RunAnnoRSCript = true ] && [ $E4SpliceCalls = false ]; then
	# Dette script tager selv fat i nødvendige Fastfiler fra MultiFQFile
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_v2FromGenotypeCalls.R $MainF $MainF/Annotation_results $SuperRunName $MultiFQFile
	# Fix For No call script, hvor multiple aminosyrer ændringer vises med ", " og ændrer til ","
	sed -i 's/,\s/,/g' $MainF/Annotation_results/ForNoCallScript_FASTQfiles_${SuperRunName}_Nuc_change_coords.txt

fi


if [ $E4SpliceCalls = true ]; then
	# Dette script tager selv fat i nødvendige Fastqfiler fra MultiFQFile
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_wE4SpliceGeneFix.R $MainF $MainF/Annotation_results $SuperRunName $MultiFQFile
	# Fix For No call script, hvor multiple aminosyrer ændringer vises med ", " og ændrer til ","
	sed -i 's/,\s/,/g' $MainF/Annotation_results/ForNoCallScript_FASTQfiles_${SuperRunName}_Nuc_change_coords.txt
fi



#  Dette script skal bruge liste over alle referencer der køres 
#  Fungerer kun hvis alle fastq filer har været alignet til samme referencer
if [ $RunSiteCovScript = true ]; then
	START=1
	END=$(awk 'END{print NR}' $RunFiles)
	for (( linenumber=$START; linenumber<=$END; linenumber++ )); do
	FastqInput=$(awk "NR==$linenumber" $RunFiles)
	FastqRunInput=$(awk "NR==$linenumber" $RunFilesR) 
	"${BashScriptF}"/Specific_site_covFromCovGenotype.sh $MainF $FastqRunInput $FastqInput $SuperRunName
	done
fi



# Kører R No_calls script 

if [ $RunNoCallsScript = true ]; then

	# Laver 1 fil med alle mulige referencer til No_call_script
	find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
	sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_latest.txt

	#if [ -d $ResultsF/$FastqRunInput/$refType ]; then
	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF $MainF/Annotation_results $SuperRunName $MultiFQFile
	# [SaveDir] [ReferenceName] [SuperRunName] [ListOfFastqFiles]

	# Samler splittede annotationresultat filer
	# Undgår første linje, da den er kolonnenavne
	rm -f $MainF/Annotation_results/AnnotationFrequency_FASTQfiles_${SuperRunName}_AllResults.txt
	cat $MainF/Annotation_results/AnnotationFrequency_FASTQfiles_${SuperRunName}*.txt | awk "NR==1 {print}" > colname_${SuperRunName}.txt
	awk FNR!=1 $MainF/Annotation_results/AnnotationFrequency_FASTQfiles_${SuperRunName}*.txt > anno_${SuperRunName}.txt 
	cat colname_${SuperRunName}.txt anno_${SuperRunName}.txt > $MainF/Annotation_results/AnnotationSummary_${SuperRunName}.txt
	rm $MainF/Annotation_results/AnnotationFrequency_FASTQfiles_${SuperRunName}*.txt
	rm colname_${SuperRunName}.txt anno_${SuperRunName}.txt
	cat $MainF/Annotation_results/FASTQfiles_${SuperRunName}_Nuc_change_coords_*.bed > $MainF/Annotation_results/VariantPositions_${SuperRunName}.bed
	rm $MainF/Annotation_results/FASTQfiles_${SuperRunName}_Nuc_change_coords_*.bed

	# Flytter filer for convenience
	cp $MainF/Annotation_results/AnnotationSummary_${SuperRunName}.txt $MainF/Results/$SuperRunName/AnnotationSummary_${SuperRunName}.txt
	cp $MainF/Annotation_results/AnnotationIndividualFiles_${SuperRunName}.txt $MainF/Results/$SuperRunName/AnnotationIndividualFiles_${SuperRunName}.txt
fi

#  Så er bunden nået. Tak fordi du læste med :)