#!/bin/bash

# HPV_T${typeTier}typing main script, that controls and calls all necessary functions

function parse_yaml {
	#  Enables yaml file being turned into objects
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

MainF=$(dirname "$(readlink -f "$0")")
eval $(parse_yaml ${MainF}/Configurations.yaml)  # Read all yaml objects into environment

# Built in operations. Must be set after parse_yaml function.
set -u  # Exit if unset variable is called
set -e  # Exit if something exits with a non-zero status

FQF=$MainF/FASTQ
# Setting the name of the run
if [ $customName = true ]; then
	TopRunName=$cName
	echo "Using custom name and saved reference info..."
else
	Date=$(date +"%H%M_%d%m%Y") 
	TopRunName=$(echo ${RunName}_${Date}) #  append date and time to name
	find $FQF/ -maxdepth 1  -name '*.fastq' > $MainF/FASTQ/FASTQfiles_${TopRunName}.txt #  Listing fastq names in a file
fi

# Create needed folders
mkdir -p $MainF/{GenotypeCalls,References/{IndexedRef,Combined_refs},Results/$TopRunName/,FASTQ,QC/{Qual,Depth,Flagstats,DuplicateMetrics,Errors}}  


# Defining locations
FQList=$MainF/FASTQ/FASTQfiles_${TopRunName}.txt
QCF=$MainF/QC; 
RefF=$MainF/References; 
ErrorF=$QCF/Errors; 
InputRefs=$RefF/InputRefs;
Rscriptfolder=$MainF/Scripts/R;
BashScriptF=$MainF/Scripts/Bash;
PythonScriptF=$MainF/Scripts/Python;
ResultsF=$MainF/Results/$TopRunName;


########### MODULES ###########
########### Primer sequence cut (BETA) ###########
if [ $cutOutsideAmplicons = true ]; then
	startN=1
	endN=$(awk 'END{print NR}' $FQList)
	# One by one, cut primer sequences
	for (( lineN=$startN; lineN<=$endN; lineN++ )); do
	FQAddr=$(awk "NR==$lineN" $FQList)
	FQName=$(basename $FQAddr .fastq)
	"${BashScriptF}"/PrimerSeqCut.sh -s "${TopRunName}" -m "$MainF" -a "$AmpliconRef" -b "$BedFileNameX"  
 	done
fi

########### FASTQ FILTRERING ###########
if [ $qualityFilt = true ]; then
	startN=1
	endN=$(awk 'END{print NR}' $FQList)
	# One by one, filter with cutadapt
	for (( lineN=$startN; lineN<=$endN; lineN++ )); do
		FQAddr=$(awk "NR==$lineN" $FQList)
		FQName=$(basename $FQAddr .fastq)
		mkdir -p $MainF/FASTQ/${TopRunName}_filtered
		# This gives filtered files the ending _filt.fastq
		cutadapt -q $QualTrim -m $MinLen -M $MaxLen "${FQAddr}" -o $FQF/${TopRunName}_filtered/"${FQName}".fastq
	done
	# Update the list of files with filt names
	find $FQF/${TopRunName}_filtered/ -maxdepth 1  -name '*.fastq' \
	> $FQF/FASTQfiles_${TopRunName}.txt    # Finding FASTQ files in FASTQ folder and create text file with them all
	FQF=$FQF/${TopRunName}_filtered # Updating varible for filtered files
	echo Done filtering fastq files
fi

########### REFERENCE INDEXING OF ALL AVAILABLE REFERENCES ###########
# Index all references in References folder. Individually puts them in 
# a subfolder in References/IndexedRef and generates variuos dict files 
# needed.
if [ $indexReferences = true ]; 
then
	for f in $RefF/*.fasta;
	do
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
		if [ ! -f ${Ref_FASTA%.fasta}.dict ]; then
			java -jar ~/picard.jar CreateSequenceDictionary \
			R=$Ref_FASTA \
			O=${Ref_FASTA%.fasta}.dict 2> $ErrorF/CreateSequenceDictionary_errors_${refType}.txt
		fi
		bwa index $Ref_FASTA
		fi
	done
fi

# Ignores Inputrefs and VirStrain subtyping. All fastq files will be aligned to the set reference in the cRef option. 
if [ $customRefForAll = true ]; 
then
	startN=1
	endN=$(awk 'END{print NR}' $FQList)
	for (( lineN=$startN; lineN<=$endN; lineN++ )); 
	do
		FQAddr=$(awk "NR==$lineN" $FQList)
		FQName=$(basename $FQAddr .fastq)
		rm -f $InputRefs/${FQName}.txt
		touch $InputRefs/${FQName}.txt
		for newref in ${cRef[@]}; 
		do
			echo $newref >> $InputRefs/${FQName}.txt
		done
	done
fi

rm -f $RefF/*fasta.fai # Can be removed? Think it secures indexing does not look on fai files and crashes 

########### GENO- AND SUBTYPING ###########
########### COVERAGE INFO GENOTYPING ###########
if [ $CovMatrixGenoTyping = true ]; then
	# Setting tier of type level that will be output
	typeTier=1
	python $PythonScriptF/autoDetectTypeFromCov.py -i $MainF/FASTQ/covMatrix.csv -o $InputRefs -r $MainF/GenotypeCalls -t $typeTier
fi


########### VIRSTRAIN GENOTYPING ###########
if [ $VirStrainGenoTyping = true ]; 
then
	# Setting tier of type level that will be output
	typeTier=1
	# Finding top 3 most possible strains, using top 1 for alignment
	VSdb=$MainF/References/$VirStrainMaindb
	startN=1
	endN=$(awk 'END{print NR}' $FQList)
	for (( lineN=$startN; lineN<=$endN; lineN++ ))
	do
		FQAddr=$(awk "NR==$lineN" $FQList)
		"${BashScriptF}"/VirStrain_genotyping.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -t ${typeTier}
		echo $lineN of $endN
	done
	# Moving summary to results folder
	cp $MainF/VirStrain_run/$TopRunName/GenotypeCalls/VirStrain_summary.txt $MainF/Results/$TopRunName/VirStrainGenotypeCalls.txt
fi


########### VIRSTRAIN SUBTYPING  ###########
# genotyped with covMatrix or manually entered in text files in References/InputRefs
if [ $VirStrainSubTyping = true ]; then
# Setting tier of type level that will be input
typeTier=1
startN=1
endN=$(awk 'END{print NR}' $FQList)
for (( lineN=$startN; lineN<=$endN; lineN++ ))
do
	FQAddr=$(awk "NR==$lineN" $FQList)
	FQName=$(basename $FQAddr .fastq)
	#  Getting correct VirStrain subtypedatabase for genotype
	outFolder=$MainF/GenotypeCalls/$TopRunName/${FQName}
	# Saving main calls to GenotypeCalls Folder
	MainCallFile=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}.txt
	#  For each genotype in inputref file, find subtype
	genoStart=1
	genoEnd=$(awk 'END{print NR}' "$MainCallFile")
	for (( genonum=$genoStart; genonum<=$genoEnd; genonum++ )); 
	do
		VirStrainSub=$(sed -n "${genonum}"p "$MainCallFile" | grep -o "^HPV[0-9]*") #  -m to only search defined line and -o to only output match
		MainCallFullName=$(sed -n "${genonum}"p "$MainCallFile") # | grep "^HPV[0-9]*"
		echo Finding $VirStrainSub subtype for ${FQName}...
		BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
		VirStrainSubDB=$MainF/References/$VirStrainSubFolders/$BaseDBName
		"${BashScriptF}"/VirStrain_subtyping.sh -f $FQAddr -s $TopRunName -m $MainF -v $VirStrainSubDB -c $MainCallFullName -t ${typeTier}
	done
	# Moving summary to results folder
	cp $MainF/VirStrain_run/$TopRunName/VirStrain_summary.txt $MainF/Results/$TopRunName/VirStrainSubtypeCalls.txt
done
fi

########### Combining references and/or saving in genotypecalls folder ###########
mkdir -p $MainF/GenotypeCalls/$TopRunName
if [ $CombineRefs = true ]; 
then
	# Define which typetier to use
	typeTier=2
	startN=1
	endN=$(awk 'END{print NR}' $FQList)
	for (( lineN=$startN; lineN<=$endN; lineN++ )); 
	do
		FQAddr=$(awk "NR==$lineN" $FQList)
		FQName=$(basename $FQAddr .fastq)
		# Finding number of references in reference file
		RefNumber=$(wc -l "$MainF"/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt | awk '{ print $1 }')
		# If more than one type in sample, combine reference files (fasta) and save info in GenotypeCalls folder
		if [ $RefNumber -gt 1 ]; then
			# Read references as array
			readarray -t RefsArr < $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt
			# Initiate array
			declare -a AllAddressList="" 
			for ref in ${RefsArr[@]}; do
				# Find ref path and append to list
				AllAddressList+=" "
				AllAddressList+=$( find $RefF -maxdepth 1  -name "${ref}*" -type f )
				#echo "${ref}" >> $MainF/GenotypeCalls/${TopRunName}/${FQName}_SplitTo.txt
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
			if [ ! -f "${Ref_FASTA%.fasta}".dict ];
			then
				java -jar ~/picard.jar CreateSequenceDictionary \
				R=$Ref_FASTA \
				O="${Ref_FASTA%.fasta}".dict 2> "$ErrorF"/CreateSequenceDictionary_errors_"${RefName}".txt
			fi
			# Create index for bwa mem
			bwa index $Ref_FASTA
			# Index with samtools faidx
			samtools faidx $Ref_FASTA
		else
			# Else if only one detected hpvtype in fastq, save reference from txt file in variable
			Ref_FASTA=$(< $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt)
		fi
	done
fi

###########
# Output: $MainF/GenotypeCalls/${TopRunName}/${FQName}.txt and _SplitTo.txt version

########### MULTIPLE TYPES AWARE ALIGNMENT & VARIANT CALLING ###########
# Input: References from $MainF/GenotypeCalls/$TopRunName/${FastQFile}.txt and fastq files 
# This will align each fastq to each of the found types one by one and then to all found types.
if [ $AlignAndVarCall = true ]; then 
# Set which type tier to align
typeTier=2
startN=1
endN=$(awk 'END{print NR}' $FQList)
for (( lineN=$startN; lineN<=$endN; lineN++ ))
do
	FQAddr=$(awk "NR==$lineN" $FQList)
	# Check if using multiple types for each fastq
	if [ $CombineRefs = true ]; 
	then
		GenoTypeCallFile=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}.txt
		# Get each type found for each fastq
		GenoTypeCallFile=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt 
		NSEQ=$(awk 'END{print NR}' $GenoTypeCallFile)
		for (( l=$startN; l<=$NSEQ; l++ )); 
		do
			Reference=$(awk "NR==$l" $GenoTypeCallFile)
	 		"${BashScriptF}"/AlignAndVariantCall.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier -b $BedFileNameX -a $AmpliconRef
	 		echo File $FQName with $Reference done with alignment and variant calling... $lineN of $endN fastq files
		done
	else
		GenoTypeCallFile=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt
	fi
	
	Reference=$(awk "NR==1" $GenoTypeCallFile)
	"${BashScriptF}"/AlignAndVariantCall.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier -b $BedFileNameX -a $AmpliconRef
	echo File $FQName with $Reference "done" with merged file alignment and variant calling...
	
done
fi

########### FOR EACH UNIQUE CHRNAME: SPLIT & VARIANTCALL AGAIN ###########
# This will split merged bamfiles into multiple bam files if they were aligned to multiple references at once.

if [ $CombineRefs = true ]; then
startN=1
endN=$(awk 'END{print NR}' $FQList)
typeTier=2
for (( lineN=$startN; lineN<=$endN; lineN++ ))
do
	FQAddr=$(awk "NR==$lineN" $FQList)
	FQName=$(basename $FQAddr .fastq)

	# Now combine maintypes to split reads into different maintypes, if there is any 
	GenoTypeCallIn=$MainF/GenotypeCalls/$TopRunName
	GenoTypeCallForFastq=$(< $GenoTypeCallIn/${FQName}_T${typeTier}.txt)

	#  Collect each split file in a folder and call variants
	StartRefs=1
	EndRefs=$(awk 'END{print NR}' $GenoTypeCallIn/${FQName}_T${typeTier}_SplitTo.txt)
	if [ $EndRefs -gt 1 ]; then 
		BamName=${FQName}_${GenoTypeCallForFastq}.sort.bam
		bamtools split -in $ResultsF/$FQName/$GenoTypeCallForFastq/ResultFiles/$BamName -reference
		for (( linen=$StartRefs; linen<=$EndRefs; linen++ )); do
			currentRef=$(awk "NR==$linen" $GenoTypeCallIn/${FQName}_T${typeTier}_SplitTo.txt)
			mkdir -p $ResultsF/$FQName/$currentRef
			SplitBam=$ResultsF/$FQName/$GenoTypeCallForFastq/ResultFiles/${FQName}_${GenoTypeCallForFastq}.sort.REF_${currentRef}.bam
			SplitBamUn=$ResultsF/$FQName/$GenoTypeCallForFastq/ResultFiles/${FQName}_${GenoTypeCallForFastq}.sort.REF_unmapped.bam
			
			# If succesful split, move bam to results folder and start variant call
			if [ -f $SplitBam ]; 
			then
				mv $SplitBam $ResultsF/$FQName/$currentRef/${FQName}_${currentRef}.bam
				#  Variant filtrering, indexing and sorting
				"${BashScriptF}"/VariantCallFromCovGenotypingSplitRefs.sh -f $FQName -s $TopRunName -m $MainF -g $GenoTypeCallIn -t $typeTier -l $currentRef -b $BedFileNameX -a $AmpliconRef 
			fi

			if [ -f $SplitBamUn ]; 
			then
				mv $SplitBamUn $ResultsF/$FQName/$currentRef/${FQName}_unmapped.bam
			fi
		done
	fi
	echo $FQName done with split calls... $lineN of $endN fastq files
done
fi
###########


###########

# Fixing name of vcf files that have had all filtered variants removed if there were no combined references with split bamfiles 
startN=1
endN=$(awk 'END{print NR}' $FQList)
typeTier=2
for (( lineN=$startN; lineN<=$endN; lineN++ ))
do
	FQAddr=$(awk "NR==$lineN" $FQList)
	FQName=$(basename $FQAddr .fastq)

	FILE1=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt
	FILE2=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}.txt

	# Checking if a fastqfile had combined references or not
	if cmp --silent -- "$FILE1" "$FILE2"; 
	then
		StartRefs=1
		EndRefs=$(awk 'END{print NR}' $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt)
		for (( linen=$StartRefs; linen<=$EndRefs; linen++ )); 
		do
			currentRef=$(awk "NR==$linen" $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt)
			# Renaming so R script can find edited vcf file
			if [ -f $ResultsF/$FQName/${currentRef}/${currentRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf ];
			then
				cp $ResultsF/$FQName/${currentRef}/${currentRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf $ResultsF/$FQName/${currentRef}/${FQName}_${currentRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf
			fi
		done
	fi
done


###########  Annotate vcfs ###########
# if [ $RunAnnoRSCript = true ] && [ $E4SpliceCalls = false ]; 
# then
# 	# Dette script tager selv fat i nødvendige Fastqfiler fra MultiFQName
# 	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_v2FromGenotypeCalls.R $MainF $MainF/Annotation_results $TopRunName $FQList
# 	# Fix For No call script, hvor multiple aminosyrer ændringer vises med ", " og ændrer til ","
# 	sed -i 's/,\s/,/g' $MainF/Annotation_results/ForNoCallScript_FASTQfiles_${TopRunName}_Nuc_change_coords.txt
# fi

########### E4 Splice calls ###########
if [ $AnnotateVariants = true ]; then
	typeTier=2
	# Script gets fastq files from MultiFQName
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_wE4SpliceGeneFix.R $MainF $MainF/Annotation_results $TopRunName $FQList $typeTier
	# Fix For No call script, hvor multiple aminosyrer ændringer vises med ", " og ændrer til ","
	sed -i 's/,\s/,/g' $MainF/Annotation_results/ForNoCallScript_${TopRunName}_Nuc_change_coords.txt
fi


########### Site coverage ###########
#  Fungerer kun hvis alle fastq filer har været alignet til samme referencer
# if [ $RunSiteCovScript = true ]; then
# 	startN=1
# 	endN=$(awk 'END{print NR}' $FQList)
# 	for (( lineN=$startN; lineN<=$endN; lineN++ )); do
# 		FQAddr=$(awk "NR==$lineN" $FQList)
# 		FQName=$(basename $FQAddr .fastq)
# 		"${BashScriptF}"/Specific_site_covFromCovGenotype.sh $MainF $FQName $TopRunName
# 	done
# fi


# ########### No_calls script ###########
if [ $Summarize = true ]; then

	# Laver 1 fil med alle mulige referencer til No_call_script
	#find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
	#sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_latest.txt

	#if [ -d $ResultsF/$FQName/$refType ]; then
	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF $MainF/Annotation_results $TopRunName $FQList
	# [SaveDir] [ReferenceName] [TopRunName] [ListOfFastqFiles]

	# Samler splittede annotationresultat filer
	cat $MainF/Annotation_results/AnnotationFrequency_${TopRunName}*.txt | awk "NR==1 {print}" > colname_${TopRunName}.txt
	awk FNR!=1 $MainF/Annotation_results/AnnotationFrequency_${TopRunName}*.txt > anno_${TopRunName}.txt  # Undgår første linje, da den er kolonnenavne
	cat colname_${TopRunName}.txt anno_${TopRunName}.txt > $MainF/Annotation_results/AnnotationSummary_${TopRunName}.txt
	rm $MainF/Annotation_results/AnnotationFrequency_${TopRunName}*.txt
	rm colname_${TopRunName}.txt anno_${TopRunName}.txt
	cat $MainF/Annotation_results/${TopRunName}_Nuc_change_coords_*.bed > $MainF/Annotation_results/VariantPositions_${TopRunName}.bed
	rm $MainF/Annotation_results/${TopRunName}_Nuc_change_coords_*.bed

	# Moving summaries to results folder
	cp $MainF/Annotation_results/AnnotationSummary_${TopRunName}.txt $MainF/Results/$TopRunName/AnnotationSummary_${TopRunName}.txt
	cp $MainF/Annotation_results/AnnotationIndividualFiles_${TopRunName}.txt $MainF/Results/$TopRunName/AnnotationIndividualFiles_${TopRunName}.txt
fi

# Cleanup
if [ $DebugMode = false ]; 
then
	if [ ${#RunName} -gt 0 ] && [ ${#MainF} -gt 0 ]; 
	then
		rm $FQList
		rm -r $FQF/${RunName}_filtered
		rm Annotation_results/*.bed Annotation_results/*.txt 
		rm -r $MainF/GenotypeCalls/${RunName}*
		rm -r $MainF/VirStrain_run/${RunName}*
	fi
fi

#  You reached the end :)