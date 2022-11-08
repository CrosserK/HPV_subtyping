#!/bin/bash

# HPV subtyping main script, that controls and calls all necessary functions

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

# Set main folder 
MainF=$(dirname "$(readlink -f "$0")")
# Read all yaml objects into environment
eval $(parse_yaml ${MainF}/Configurations.yaml)  

# Built in operations. Must be set after parse_yaml function.
set -u  # Exit if unset variable is called
set -e  # Exit if something exits with a non-zero status

# Set input folder
FQF=$MainF/FASTQ

############################################
# Check configurations
if [ $CovMatrixGenoTyping = true ] && [ $VirStrainGenoTyping = true ]; then
	echo WARNING! Both genotyping options are turned on in Configurations.yaml. Choose only one and start again. 
	exit 1
else
	echo Starting FindHPVMutations program...
fi

if [ $CovMatrixGenoTyping = true ]; then
	if [ ! -f $FQF/covMatrix.csv ]; then
		echo WARNING! No covMatrix.csv found in he FASTQ folder. Exiting...
		exit 1
	fi
fi

if [ $startPreviousRun = true ]; then
	read -p "WARNING. Continuing previous run. Are you sure you wanted this? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
fi
###############################################

# Setting the name of the run (TopRunName)
if [ $startPreviousRun = false ]; then
	#  append date and time to given name
	Date=$(date +"%H%M_%d%m%Y") 
	TopRunName=$(echo ${RunName}_${Date}) 
	#  Listing fastq names in a file
	find $FQF/ -maxdepth 1  -name '*.fastq' > $FQF/FASTQfiles_${TopRunName}.txt
else
	# This assumes the FQList already exist, because the only reason to use startPreviousRun, is to continue a previous run.
	TopRunName=$(ls -Art $FQF | grep "FASTQfiles_" | tail -n 1 | sed 's/FASTQfiles_//g' | sed 's/.txt//g') # This gets the lafinding run name
	echo "Starting previous run and using saved fastq info..."
fi

# Create needed folders
mkdir -p $MainF/{References/{IndexedRef,Combined_refs},Results/$TopRunName/,FASTQ,QC/{Flagstats,Logs}}  

# Defining locations
FQList=$FQF/FASTQfiles_${TopRunName}.txt
QCF=$MainF/QC; 
RefF=$MainF/References; 
LogF=$QCF/Logs; 
InputRefs=$RefF/InputRefs;
Rscriptfolder=$MainF/Scripts/R;
BashScriptF=$MainF/Scripts/Bash;
PythonScriptF=$MainF/Scripts/Python;
ResultsF=$MainF/Results/$TopRunName;

# Set number of loops from FQList. Paths will be gathered from FQlist file
FQListLen=$(awk 'END{print NR}' $FQList)

########### MODULES ###########

########### FASTQ FILTRERING ###########
if [ $qualityFilt = true ]; then
	# One by one, filter with cutadapt
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ )); 
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)
		mkdir -p $FQF/${TopRunName}_filtered
		# This gives filtered files the ending '_filt.fastq'
		cutadapt -q $QualTrim -m $MinLen -M $MaxLen "${FQAddr}" \
		-o $FQF/${TopRunName}_filtered/"${FQName}".fastq
	done
	# Update the list of files with '_filt' names. 
	find $FQF/${TopRunName}_filtered/ -maxdepth 1 -name '*.fastq' \
	> $FQF/FASTQfiles_${TopRunName}.txt  
	FQF=$FQF/${TopRunName}_filtered
	echo Done filtering fastq files
fi

########### REFERENCE INDEXING OF ALL AVAILABLE REFERENCES ###########
# Index all references in References folder. Individually puts them in 
# a subfolder in References/IndexedRef and generates variuos index files.
if [ $indexReferences = true ]; 
then
	for f in $RefF/*.fasta;
	do
		RefName=${f##*/}
		mkdir -p $RefF/IndexedRef/${RefName%.fasta}

		# Samtools indexing
		if [ -f $RefF/IndexedRef/${RefName%.fasta}/${RefName} ]; then
			echo $RefName already in IndexedRef folder
		else
			cp $f $RefF/IndexedRef/${RefName%.fasta}/$RefName
			Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName
			samtools faidx $Ref_FASTA
		# GATK indexing haplotypecaller
		if [ ! -f ${Ref_FASTA%.fasta}.dict ]; then
			java -jar ~/picard.jar CreateSequenceDictionary \
			R=$Ref_FASTA \
			O=${Ref_FASTA%.fasta}.dict 2> $LogF/CreateSequenceDictionary_log_${refType}.txt
		fi
		# BWA indexing
		bwa index $Ref_FASTA
		fi
	done
fi

if [ $sameRefForAll = true ]; 
# Ignores Inputrefs and VirStrain subtyping. All fastq files will be aligned to the set reference in the sRef option. 
then
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ )); 
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)
		rm -f $InputRefs/${FQName}.txt
		# For each listed ref append to new line in input ref for fastq
		for newref in ${sRef[@]}; 
		do
			echo $newref >> $InputRefs/${FQName}.txt
		done
	done
fi

########### COVERAGE INFO GENOTYPING ###########
# Find the genotypes for each fastq from a coverage matrix from Ion Torrent sequencing. 
if [ $CovMatrixGenoTyping = true ]; then
	# Setting tier of type level that will be output (etc. genotype = 1, subtype = 2)
	typeTier=1
	python $PythonScriptF/autoDetectTypeFromCov.py -i $MainF/FASTQ/covMatrix.csv -o $ResultsF -r $RefF -t $typeTier
fi


########### VIRSTRAIN GENOTYPING ###########
if [ $VirStrainGenoTyping = true ]; 
then
	# Setting tier of type level that will be output (etc. genotype = 1, subtype = 2)
	typeTier=1
	# Main data base fetched from yaml file
	VSdb=$MainF/References/$VirStrainMaindb
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ ))
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)
		echo Genotyping $FQName
		"${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -t ${typeTier} -g "$VirSGenoTypingTop"
		echo "$FQNum of $FQListLen done"
	done
fi


########### VIRSTRAIN SUBTYPING  ###########
#  Getting correct VirStrain subtypedatabase for genotype
#  Saving main calls to GenotypeCalls Folder
if [ $VirStrainSubTyping = true ]; then
# Setting tier of type rank that will be input (if it is best match or 2nd best match etc.)
typeRank=1
for (( FQNum=1; FQNum<=$FQListLen; FQNum++ ))
do
	FQAddr=$(awk "NR==$FQNum" $FQList)
	FQName=$(basename $FQAddr .fastq)

	typeTier=1
	MainCallFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_SplitTo.txt
	#  Setting typeTier to output by incrementing input typeTier to output typeTier
	let "typeTier=typeTier+1"
	#  For each genotype in inputref file, find subtype
	if [ -f $MainCallFile ]; # Skips fastq if no genotype were found. 
	then
		:
	else
		continue
	fi

	genoEnd=$(awk 'END{print NR}' "$MainCallFile")
	for (( genonum=1; genonum<=$genoEnd; genonum++ )); 
	do
		VirStrainSub=$(sed -n "${genonum}"p "$MainCallFile" | grep -o "^HPV[0-9]*") # -o to only output match
		prevCall=$(sed -n "${genonum}"p "$MainCallFile" | grep "^HPV[0-9]*" | awk '{print $1}')
		# varCheck() {
		# 	# Checks if variable is defined, else skips loop 
		# 	if [ -n "${prevCall-}" ]
		# 	then
		# 		echo passing
		# 		:
		# 	elif [ "${prevCall+defined}" = defined ]
		# 	then
		# 		continue
		# 	else
		# 		continue
		# 	fi
		# }
		#varCheck
		echo Finding $VirStrainSub subtype for ${FQName}...
		BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
		VSdb=$MainF/References/$VirStrainSubFolders/$BaseDBName # VirStrainSubFolders comes from Configurations.yaml file
		echo "${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -c $prevCall -t $typeTier -g "$VirSSubTypingTop" -r "$genonum"
		"${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -c $prevCall -t $typeTier -g "$VirSSubTypingTop" -r "$genonum"
		echo "$FQNum of $FQListLen done"
	done
done
fi

# ########### Combining references and/or saving in genotypecalls folder ###########
if [ $CombineRefs = true ]; 
then
	# Setting tier of type level that will be used (etc. genotype = 1, subtype = 2)
	typeTier=2
	# Setting tier of type rank that will be input (if it is best match (1) or 2nd best match (2) etc.)
	typeRank=1
	# removing unwanted fai files
	rm -f $RefF/*.fasta.fai
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ )); 
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)
		# Finding number of references in reference file
		# If more than one type in sample, combine reference files (fasta) and save info in GenotypeCalls folder
		if [ $VirSSubTypingTop -gt 1 ]; then
			# Read references as array
			readarray -t RefsArr < $ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}_SplitTo.txt
			# Initiate array
			declare -a AllAddressList="" 
			for ref in ${RefsArr[@]}; do
				# Find ref path and append to list
				echo $ref
				AllAddressList+=" "
				AllAddressList+=$( find $RefF -maxdepth 1  -name "${ref}*" -type f )
			done
			# Make new fasta name for combined ref
			newRefName=$(printf "_%s" "${RefsArr[@]}")
			newRefName=${newRefName:1} # Removes wrong prefix of "_" from above method
			newRefNameFasta=${newRefName}.fasta
			# Collect the fastas to one file and index, if it does not exist
			if [ ! -f $RefF/Combined_refs/$newRefNameFasta ];
			then
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
				O="${Ref_FASTA%.fasta}".dict 2> "$LogF"/CreateSequenceDictionary_log_"${RefName}".txt
				# Create index for bwa mem
				bwa index $Ref_FASTA
				# Index with samtools faidx
				samtools faidx $Ref_FASTA
			fi
		fi
	done
fi


########### ALIGNMENT & VARIANT CALLING (INCL. MERGED REF) ###########
# This will align each fastq to each of the found types one by one and then to all found types.
if [ $AlignAndVariantCall = true ]; 
then 
	#  Set which type tier to align
	typeTier=2
	typeRank=1
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ ));
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)
		#  Check if using combined ref
		if [ $CombineRefs = true ]; 
		#  Align to combined reference
		then
			#  Get each type found for each fastq
			GenoTypeCallFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}.txt 
			Reference=$(awk "NR==1" $GenoTypeCallFile)
			"${BashScriptF}"/Alignment.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier
			echo File $FQName with $Reference "done" with MERGED file alignment... 
			
			#  Get merged ref name
			MergedRef=$(< $ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}.txt)
			EndRefs=$(awk 'END{print NR}' $ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}_SplitTo.txt)
			#  Split the combined reference alignments into individual bams and variant call
			if [ $EndRefs -gt 1 ];
			then 

				#  Split bam
				BamName=${FQName}_${MergedRef}.sort.bam
				bamtools split -in $ResultsF/$FQName/$MergedRef/ResultFiles/$BamName -reference

				#  Collect each split file in a folder and call variants
				for (( RefNum=1; RefNum<=$EndRefs; RefNum++ )); 
				do
					
					currentRef=$(awk "NR==$RefNum" $ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}_SplitTo.txt)
					mkdir -p $ResultsF/$FQName/$currentRef/ResultFiles
					SplitBam=$ResultsF/$FQName/$MergedRef/ResultFiles/${FQName}_${MergedRef}.sort.REF_${currentRef}.bam
					SplitBamUn=$ResultsF/$FQName/$MergedRef/ResultFiles/${FQName}_${MergedRef}.sort.REF_unmapped.bam
					
					#  If succesful split, move bam to results folder and start variant call
					if [ -f $SplitBam ]; 
					then
						# Rename from what bamtools split called them
						mv $SplitBam $ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_${currentRef}.sort.bam
						SplitBam=$ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_${currentRef}.sort.bam

						#  Variant filtrering, indexing and sorting
						echo -f $FQName -s $TopRunName -m $MainF -b $SplitBam -t $typeTier -r $typeRank -l $currentRef
						"${BashScriptF}"/VariantCall.sh -f $FQName -s $TopRunName -m $MainF -r $currentRef -b $SplitBam -t $typeTier -a "ampliconRefPlaceholder"
						echo File $FQName with $Reference "done" with split file variant call...
					fi

					if [ -f $SplitBamUn ]; 
					then
						# Rename from what bamtools split called them
						mv $SplitBamUn $ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_unmapped.sort.bam
						SplitBamUn=$ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_unmapped.sort.bam
					fi
				done
			else
				# Align without splits
				GenoTypeCallFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}_SplitTo.txt
				Reference=$(awk "NR==1" $GenoTypeCallFile)
				"${BashScriptF}"/Alignment.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier
				"${BashScriptF}"/Alignment.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier -b "bedfileReplacement"
				
				# Variant call
				workD=$ResultsF/$FQName
				currentF=$workD/$Reference
				BamFile=$currentF/${Reference}.bam
				BamFile="${BamFile%bam}"sort.bam
				"${BashScriptF}"/VariantCall.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -b $BamFile -t $typeTier -a "ampliconRefPlaceholder"
				echo File $FQName with $Reference "done" with alignment and variant calling... $FQNum of $FQListLen fastq files
			fi
		echo $FQName "done" with split calls... $FQNum of $FQListLen fastq files
		else
			# Align without splits
	 		GenoTypeCallFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}_SplitTo.txt
			Reference=$(awk "NR==1" $GenoTypeCallFile)
	 		"${BashScriptF}"/Alignment.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier -b "bedfileReplacement"
			
			# Variant call
			workD=$ResultsF/$FQName
			currentF=$workD/$Reference
			BamFile=$currentF/${Reference}.bam
			BamFile="${BamFile%bam}"sort.bam
			"${BashScriptF}"/VariantCall.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -b $BamFile -t $typeTier -a "ampliconRefPlaceholder"
	 		echo File $FQName with $Reference "done" with alignment and variant calling... $FQNum of $FQListLen fastq files
		fi
	done
fi



###########

# Fixing name of vcf files that have had all filtered variants removed if there were no combined references with split bamfiles 
typeTier=2
for (( FQNum=1; FQNum<=$FQListLen; FQNum++ ))
do
	FQAddr=$(awk "NR==$FQNum" $FQList)
	FQName=$(basename $FQAddr .fastq)

	FILE1=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt
	FILE2=$MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}.txt

	# Checking if a fastqfile had combined references or not
	if cmp --silent -- "$FILE1" "$FILE2"; 
	then
		StartRefs=1
		EndRefs=$(awk 'END{print NR}' $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt)
		for (( FQNum=$StartRefs; FQNum<=$EndRefs; FQNum++ )); 
		do
			currentRef=$(awk "NR==$FQNum" $MainF/GenotypeCalls/$TopRunName/${FQName}/${FQName}_T${typeTier}_SplitTo.txt)
			# Renaming so R script can find edited vcf file
			if [ -f $ResultsF/$FQName/${currentRef}/${currentRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf ];
			then
				cp $ResultsF/$FQName/${currentRef}/${currentRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf $ResultsF/$FQName/${currentRef}/${FQName}_${currentRef}.sort.readGroupFix_filtered_FiltEx_headerfix.vcf
			fi
		done
	fi
done


########### Annotate variants including E4 Splice calls ###########
if [ $AnnotateVariants = true ]; then
	typeTier=2
	# Script gets fastq files from MultiFQName
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_wE4SpliceGeneFix.R $MainF $ResultsF $TopRunName $FQList $typeTier
	# Fix For No call script, hvor multiple aminosyrer ændringer vises med ", " og ændrer til ","
	sed -i 's/,\s/,/g' $ResultsF/ForNoCallScript_${TopRunName}_Nuc_change_coords.txt
fi


# ########### No_calls script ###########
if [ $Summarize = true ]; then

	# Laver 1 fil med alle mulige referencer til No_call_script
	#find $RefF/ -maxdepth 1 -name '*.fasta' | sed 's/^.*\(References.*fasta\).*$/\1/' | \
	#sed 's/.fasta//g' | sed 's/References\///g' > $RefdF/RefSubtyper_latest.txt

	#if [ -d $ResultsF/$FQName/$refType ]; then
	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF $ResultsF $TopRunName $FQList
	# [SaveDir] [ReferenceName] [TopRunName] [ListOfFastqFiles]

	# Samler splittede annotationresultat filer
	if [ -f $ResultsF/AnnotationFrequency_${TopRunName}*.txt ]; then
		cat $ResultsF/AnnotationFrequency_${TopRunName}*.txt | awk "NR==1 {print}" > colname_${TopRunName}.txt
		awk FNR!=1 $ResultsF/AnnotationFrequency_${TopRunName}*.txt > anno_${TopRunName}.txt  # Undgår første linje, da den er kolonnenavne
		cat colname_${TopRunName}.txt anno_${TopRunName}.txt > $ResultsF/AnnotationSummary_${TopRunName}.txt
		rm $ResultsF/AnnotationFrequency_${TopRunName}*.txt
		rm colname_${TopRunName}.txt anno_${TopRunName}.txt
		cat $ResultsF/${TopRunName}_Nuc_change_coords_*.bed > $ResultsF/VariantPositions_${TopRunName}.bed
		rm $ResultsF/${TopRunName}_Nuc_change_coords_*.bed
	else
		rm $ResultsF/ForNoCallScript_*.txt
		echo "No variants found in this run" > $ResultsF/AnnotationSummary_${TopRunName}.txt
	fi
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