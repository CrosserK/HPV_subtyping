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
mkdir -p $MainF/{References/{IndexedRef,Combined_refs},Results/$TopRunName/,FASTQ,QC/{Flagstats,Logs,MarkDup,Indexing}}  

# Defining locations
FQList=$FQF/FASTQfiles_${TopRunName}.txt
QCF=$MainF/QC; 
RefF=$MainF/References; 
LogF=$QCF/Logs;
IdxF=$QCF/Indexing; 
InputRefs=$RefF/InputRefs;
Rscriptfolder=$MainF/Scripts/R;
BashScriptF=$MainF/Scripts/Bash;
PythonScriptF=$MainF/Scripts/Python;
ResultsF=$MainF/Results/$TopRunName;

# Clear previous logs
find $QCF -type f -delete

echo [$(date +"%d-%m-%Y %H:%M:%S")] "Starting FindHPVMutations program..." >> $LogF/${TopRunName}.txt
# Set number of loops from FQList. Paths will be gathered from FQlist file
FQListLen=$(awk 'END{print NR}' $FQList)

########### MODULES ###########

########### FASTQ FILTRERING ###########
if [ $qualityFilt = true ]; then
	# One by one, filter with cutadapt
	echo [$(date +"%d-%m-%Y %H:%M:%S")] "Filtering fastq files..." >> $LogF/${TopRunName}.txt
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
	echo [$(date +"%d-%m-%Y %H:%M:%S")] Done filtering fastq files >> $LogF/${TopRunName}.txt
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
			echo [$(date +"%d-%m-%Y %H:%M:%S")] $RefName already in IndexedRef folder >> $LogF/${TopRunName}.txt
		else
			cp $f $RefF/IndexedRef/${RefName%.fasta}/$RefName
			Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName
			samtools faidx $Ref_FASTA
		# GATK indexing haplotypecaller
		if [ ! -f ${Ref_FASTA%.fasta}.dict ]; then
			java -jar ~/picard.jar CreateSequenceDictionary \
			R=$Ref_FASTA \
			O=${Ref_FASTA%.fasta}.dict 2> $IdxF/CreateSequenceDictionary_log_${refType}.txt
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
	python $PythonScriptF/autoDetectTypeFromCov.py -i $MainF/FASTQ/covMatrix.csv -o $ResultsF -r $RefF -f $FQF -t $typeTier -g $HouseGeneEndRow -l $LogF/${TopRunName}.txt
	# Integration detection by deletion detection
	python $PythonScriptF/detectIntegration.py -i $MainF/FASTQ/covMatrix.csv -o $ResultsF -r $RefF -g $HouseGeneEndRow -a $MainF/FASTQ/amplPos.bed
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
		echo [$(date +"%d-%m-%Y %H:%M:%S")] Genotyping $FQName >> $LogF/${TopRunName}.txt
		"${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -t ${typeTier} -g "$VirSGenoTypingTop"
		echo "$FQNum of $FQListLen done"
		echo [$(date +"%d-%m-%Y %H:%M:%S")] "$FQNum of $FQListLen done" >> $LogF/${TopRunName}.txt
	done
fi



# Combined ref split before subtyping
# ########### Combining references and/or saving in genotypecalls folder ###########
if [ $CombineRefs = true ] && [ $AlignAndVariantCall = true ]; 
then


	# Combine fasta refs. Removing unwanted fai files first.
	rm -f $RefF/*.fasta.fai
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ )); 
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)

		# Setting tier of type level that will be used as input (etc. genotype = 1, subtype = 2)
		typeTier=1
		# Setting tier of type rank that will be input (if it is best match = 1, or 2nd best match = 2, etc.)
		typeRank=1
		# Finding number of references in reference file
		# If more than one type in sample, combine reference files (fasta) and save info in GenotypeCalls folder
		RefFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_SplitTo.txt
		CombRefFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}.txt
		# Checking that file exists, else going to next iteration
		if [ ! -f $RefFile ]; 
		then
			echo [$(date +"%d-%m-%Y %H:%M:%S")] "No genotype file for $FQName" >> $LogF/${TopRunName}.txt
			continue
		fi

		RefFileLen=$(awk 'END{print NR}' $RefFile)
		# If more than one reference found, combine fasta and align to it, then bamsplit and variant call
		if [ $RefFileLen -gt 1 ]; 
		then

			#### CREATE COMBINED REFERENCE
			# Read references as array
			readarray -t RefsArr < $RefFile
			# Initiate array
			declare -a AllAddressList="" 
			for ref in ${RefsArr[@]}; do
				# Find ref path and append to list
				AllAddressList+=" "
				AllAddressList+=$( find $RefF -maxdepth 1  -name "${ref}*" -type f )
			done

			# Make new fasta name for combined ref
			newRefName=$(printf "_%s" "${RefsArr[@]}")
			newRefName=${newRefName:1} # Removes wrong prefix of "_" from above method
			newRefNameFasta=${newRefName}.fasta

			# Collect the fastas to one file and index, if it does not exist
			if [ ! -f $RefF/Combined_refs/${newRefName}.dict ];
			then
				cat $AllAddressList > $RefF/Combined_refs/$newRefNameFasta
				# Index and save in IndexedRefs folder
				NewRef=$RefF/Combined_refs/$newRefNameFasta
				RefName=${NewRef##*/}
				mkdir -p $RefF/IndexedRef/${RefName%.fasta}
				cp $NewRef $RefF/IndexedRef/${RefName%.fasta}/$RefName
				Ref_FASTA=$RefF/IndexedRef/${RefName%.fasta}/$RefName
				# Avoiding if there is only one ref and it is already indexed
				if [ ! -f "${Ref_FASTA%.fasta}".dict ]; then

					# Create sequence dictionary for gatk haplotypecaller
					java -jar ~/picard.jar CreateSequenceDictionary \
					R=$Ref_FASTA \
					O="${Ref_FASTA%.fasta}".dict 2> "$IdxF"/CreateSequenceDictionary_log_"${RefName}".txt
					# Create index for bwa mem
					bwa index $Ref_FASTA
					# Index with samtools faidx

					samtools faidx $Ref_FASTA
				fi
			fi

			####### COMBINED REFERENCE CREATED

			#  Align to combined reference, split and subtype
			Reference=$(awk "NR==1" $CombRefFile)
			"${BashScriptF}"/Alignment.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier
			echo File $FQName with $Reference "done" with MERGED file alignment... 
			echo [$(date +"%d-%m-%Y %H:%M:%S")] "File $FQName with $Reference "done" with MERGED file alignment... " >> $LogF/${TopRunName}.txt

			#  Get merged ref name
			MergedRef=$(< $CombRefFile)

			#  Split the combined reference alignments into individual bams and variant call
			BamName=${FQName}_${MergedRef}.sort.bam
			bamtools split -in $ResultsF/$FQName/$MergedRef/ResultFiles/$BamName -reference

			#  Collect each split file in a folder and call variants
			for (( RefNum=1; RefNum<=$RefFileLen; RefNum++ )); 
			do
				
				currentRef=$(awk "NR==$RefNum" $RefFile)
				mkdir -p $ResultsF/$FQName/$currentRef/ResultFiles
				SplitBam=$ResultsF/$FQName/$MergedRef/ResultFiles/${FQName}_${MergedRef}.sort.REF_${currentRef}.bam
				SplitBamUn=$ResultsF/$FQName/$MergedRef/ResultFiles/${FQName}_${MergedRef}.sort.REF_unmapped.bam
				
				#  If succesful split, move bam to results folder and start variant call
				if [ -f $SplitBam ]; 
				then
					# Rename from what bamtools split called them
					mv $SplitBam $ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_${currentRef}.sort.bam
					SplitBam=$ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_${currentRef}.sort.bam

					# SUBTYPE THE SPLIT BAM. For each genotype in inputref file, find subtype
					if [ $VirStrainSubTyping = true ]; then
						# Setting tier of type rank that will be input (if it is best match or 2nd best match etc.)
						typeRank=1
						typeTier=1
						#  Setting typeTier to output by incrementing input typeTier to output typeTier
						let "typeTier=typeTier+1"

						VirStrainSub=$(sed -n "${RefNum}"p "$RefFile" | grep -o "^HPV[0-9]*") # -o to only output match
						prevCall=$(sed -n "${RefNum}"p "$RefFile" | grep "^HPV[0-9]*" | awk '{print $1}')
						# Avoid edge cases where HPV is not named HPV[0-9]*, but like HPV-[A-Z]. These are not integrated into the pipeline yet (such as HPV-mSK008_MH777156_2)
						if [ $VirStrainSub = HPV ]; then
							continue
						fi
						echo Finding $VirStrainSub subtype rank $RefNum for ${FQName}...
						echo [$(date +"%d-%m-%Y %H:%M:%S")] "Finding $VirStrainSub subtype rank $RefNum for ${FQName}..." >> $LogF/${TopRunName}.txt
						BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
						VSdb=$MainF/References/$VirStrainSubFolders/$BaseDBName # VirStrainSubFolders comes from Configurations.yaml file
						"${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -c $prevCall -t $typeTier -g "$VirSSubTypingTop" -r "$RefNum"
						if [ $RefFileLen = $RefNum ]; then
							echo "$FQNum of $FQListLen "done""
							echo [$(date +"%d-%m-%Y %H:%M:%S")] "$FQNum of $FQListLen "done"" >> $LogF/${TopRunName}.txt
						fi
					fi 

					#  Variant filtrering, indexing and sorting
					"${BashScriptF}"/VariantCall.sh -f $FQName -s $TopRunName -m $MainF -r $currentRef -b $SplitBam -t $typeTier -a "ampliconRefPlaceholder"
					echo File $FQName with $Reference "done" with split file variant call...
					echo [$(date +"%d-%m-%Y %H:%M:%S")] "File $FQName with $Reference "done" with split file variant call..." >> $LogF/${TopRunName}.txt

				fi

				if [ -f $SplitBamUn ]; 
				then
					# Rename from what bamtools split called them
					mv $SplitBamUn $ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_unmapped.sort.bam
					SplitBamUn=$ResultsF/$FQName/$currentRef/ResultFiles/${FQName}_BAMSplit_unmapped.sort.bam
				fi
			done

		echo $FQName "done" with split calls... $FQNum of $FQListLen fastq files

		echo [$(date +"%d-%m-%Y %H:%M:%S")] "$FQName "done" with split calls... $FQNum of $FQListLen fastq files" >> $LogF/${TopRunName}.txt
		else
		# If only one reference subtype, align and variant call without splits

			if [ $VirStrainSubTyping = true ]; then
				# Setting tier of type rank that will be input (if it is best match or 2nd best match etc.)
				typeRank=1
				typeTier=1
				#  Setting typeTier to output by incrementing input typeTier to output typeTier
				let "typeTier=typeTier+1"

				VirStrainSub=$(sed -n 1p "$RefFile" | grep -o "^HPV[0-9]*") # -o to only output match
				prevCall=$(sed -n 1p "$RefFile" | grep "^HPV[0-9]*" | awk '{print $1}')
				# Avoid edge cases where HPV is not named HPV[0-9]*, but like HPV-[A-Z]. These are not integrated into the pipeline yet (such as HPV-mSK008_MH777156_2)
				if [ $VirStrainSub = HPV ]; then
					continue
				fi
				echo Finding $VirStrainSub subtype for ${FQName}...
				echo [$(date +"%d-%m-%Y %H:%M:%S")] "Finding $VirStrainSub subtype for ${FQName}..." >> $LogF/${TopRunName}.txt
				BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
				VSdb=$MainF/References/$VirStrainSubFolders/$BaseDBName # VirStrainSubFolders comes from Configurations.yaml file
				"${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -c $prevCall -t $typeTier -g "$VirSSubTypingTop" -r 1
				echo "$FQNum of $FQListLen "done""
				echo [$(date +"%d-%m-%Y %H:%M:%S")] "$FQNum of $FQListLen "done"" >> $LogF/${TopRunName}.txt
			fi

			# Align without splits
			typeTier=2
			typeRank=1
			RefFile=$ResultsF/${FQName}/TypeCalls/${FQName}_T${typeTier}_R${typeRank}_SplitTo.txt
			Reference=$(awk "NR==1" $RefFile)
			"${BashScriptF}"/Alignment.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -t $typeTier -b "bedfileReplacement"
			
			# Variant call
			workD=$ResultsF/$FQName
			currentF=$workD/$Reference
			BamFile=$currentF/${Reference}.bam
			BamFile="${BamFile%bam}"sort.bam
			"${BashScriptF}"/VariantCall.sh -f $FQAddr -s $TopRunName -m $MainF -r $Reference -b $BamFile -t $typeTier -a "ampliconRefPlaceholder"
			echo File $FQName with $Reference "done" with alignment and variant calling... $FQNum of $FQListLen fastq files
			echo [$(date +"%d-%m-%Y %H:%M:%S")] "File $FQName with $Reference "done" with alignment and variant calling... $FQNum of $FQListLen fastq files" >> $LogF/${TopRunName}.txt

		fi

	done
fi



########### VIRSTRAIN SUBTYPING  ###########
#  Getting correct VirStrain subtypedatabase for genotype
#  Saving main calls to GenotypeCalls Folder
if [ $VirStrainSubTyping = true ] && [ $CombineRefs = false ];	 then
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
		# Avoid edge cases where HPV is not named HPV[0-9]*, but like HPV-[A-Z]. These are not integrated into the pipeline yet (such as HPV-mSK008_MH777156_2)
		if [ $VirStrainSub = HPV ]; then
			continue
		fi
		echo Finding $VirStrainSub subtype rank $genonum for ${FQName}...
		echo [$(date +"%d-%m-%Y %H:%M:%S")] "Finding $VirStrainSub subtype rank $genonum for ${FQName}..." >> $LogF/${TopRunName}.txt
		BaseDBName=$(echo ${VirStrainSub}_VirStrainDB)
		VSdb=$MainF/References/$VirStrainSubFolders/$BaseDBName # VirStrainSubFolders comes from Configurations.yaml file
		"${BashScriptF}"/VirStrain_typeCall.sh -f $FQAddr -s $TopRunName -m $MainF -v $VSdb -c $prevCall -t $typeTier -g "$VirSSubTypingTop" -r "$genonum"
		if [ $genonum = $genoEnd ]; then
			echo "$FQNum of $FQListLen "done""
			echo [$(date +"%d-%m-%Y %H:%M:%S")] "$FQNum of $FQListLen "done"" >> $LogF/${TopRunName}.txt
		fi
	done
done
fi





########### ALIGNMENT & VARIANT CALLING (INCL. MERGED REF) ###########
# This will align each fastq to each of the found types one by one and then to all found types.
if [ $CombineRefs = false ] && [ $AlignAndVariantCall = true ]; 
then 
	#  Set which type tier to combine
	typeTier=1
	typeRank=1
	for (( FQNum=1; FQNum<=$FQListLen; FQNum++ ));
	do
		FQAddr=$(awk "NR==$FQNum" $FQList)
		FQName=$(basename $FQAddr .fastq)	

		# Align without splits
		typeTier=2
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
		echo [$(date +"%d-%m-%Y %H:%M:%S")] "File $FQName with $Reference "done" with alignment and variant calling... $FQNum of $FQListLen fastq files" >> $LogF/${TopRunName}.txt

	done
fi



########### Annotate variants including E4 Splice calls ###########
if [ $AnnotateVariants = true ]; then
	typeTier=2
	# Script gets fastq files from MultiFQName
	echo Annotating variants... 
	echo [$(date +"%d-%m-%Y %H:%M:%S")] Annotating variants... >> $LogF/${TopRunName}.txt
	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash_wE4SpliceGeneFix.R $MainF $ResultsF $TopRunName $FQList $typeTier
	# Fix For No call script, hvor multiple aminosyrer ændringer vises med ", " og ændrer til ","
	sed -i 's/,\s/,/g' $ResultsF/ForSummaryScript_${TopRunName}_Nuc_change_coords.txt
fi


# ########### Summarizing variant results ###########
if [ $Summarize = true ]; then

	cat $ResultsF/${TopRunName}_Nuc_change_coords_*.bed > $ResultsF/VariantPositions_${TopRunName}.bed
	echo Summarizing variant results...
	echo [$(date +"%d-%m-%Y %H:%M:%S")] Summarizing variant results... >> $LogF/${TopRunName}.txt
	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF $ResultsF $TopRunName $FQList
	# [SaveDir] [ReferenceName] [TopRunName] [ListOfFastqFiles]

	# Samler splittede annotationresultat filer
	if ls $ResultsF/AnnotationFrequency_${TopRunName}*.txt 1> /dev/null 2>&1; then
		# Collecting summaries from each reference
		cat $ResultsF/AnnotationFrequency_${TopRunName}*.txt | awk "NR==1 {print}" > colname_${TopRunName}.txt
		awk FNR!=1 $ResultsF/AnnotationFrequency_${TopRunName}*.txt > anno_${TopRunName}.txt  # Undgår første linje, da den er kolonnenavne
		cat colname_${TopRunName}.txt anno_${TopRunName}.txt > $ResultsF/AnnotationSummary_${TopRunName}.txt
		rm $ResultsF/AnnotationFrequency_${TopRunName}*.txt
		rm colname_${TopRunName}.txt anno_${TopRunName}.txt
		rm $ResultsF/${TopRunName}_Nuc_change_coords_*.bed
	else
		rm $ResultsF/ForSummaryScript_*.txt
		echo "No variants found in this run" > $ResultsF/AnnotationSummary_${TopRunName}.txt
		echo [$(date +"%d-%m-%Y %H:%M:%S")] "No variants found in this run. Low coverage?" >> $LogF/${TopRunName}.txt
	fi
	echo [$(date +"%d-%m-%Y %H:%M:%S")] "Summaries done" >> $LogF/${TopRunName}.txt
fi

# Cleanup
if [ $DebugMode = false ]; 
then
	if [ ${#RunName} -gt 0 ] && [ ${#MainF} -gt 0 ]; 
	then
		rm $FQList
		rm -f Annotation_results/*.bed Annotation_results/*.txt
		# Removing filtered fastq files if filtering were done
		if [[ "$FQF" == *"_filtered" ]]; then
			rm -r $FQF
		fi
	fi
fi

echo [$(date +"%d-%m-%Y %H:%M:%S")] "FindHPVMutations done." >> $LogF/${TopRunName}.txt
#  You reached the end :)