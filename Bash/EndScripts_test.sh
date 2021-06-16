
###### TEST #######
########## Define Folders ##########
SuperRunName=Exome_50_120_ampliconcalls_PaVE
MainF=/home/pato/Skrivebord/HPV16_projekt
refType=K02718.1_PaVE
SuperRunFolder=$SuperRunName
RefListFile=RefSubtyper_allFastq
####################################

RunAnnoRSCript=true
RunSiteCovSCript=true
RunNoCallsSCript=true

Rscriptfolder=$MainF/Scripts/R
MultiFQFile=$(echo FASTQfiles_$SuperRunFolder)

if [ $RunAnnoRSCript = true ]; then

	Rscriptfolder=$MainF/Scripts/R

	Rscript $Rscriptfolder/Annotate_vcf_multiple_for_bash.R $MainF/Annotation_results/ $refType $SuperRunName $MultiFQFile 

fi

# Dette script skal bruge liste over alle referencer der køres

if [ $RunSiteCovSCript = true ]; then

FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt) 
FastqRunList=$(< $MainF/FASTQfiles_${SuperRunName}_runnames.txt)

START=1
END=$(awk 'END{print NR}' $MainF/FASTQfiles_${SuperRunName}.txt)

for (( linenumber=$START; linenumber<=$END; linenumber++ ))
do
	FastqInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}.txt)
	FastqRunInput=$(awk "NR==$linenumber" $MainF/FASTQfiles_${SuperRunName}_runnames.txt) 
	Specific_site_cov.sh $FastqRunInput $FastqInput $SuperRunFolder $RefListFile
	echo $linenumber of $END
done
fi


# Kører R No_calls script 

if [ $RunNoCallsSCript = true ]; then

	Rscript $Rscriptfolder/No_calls_on_vcf.R $MainF/Annotation_results/ $refType $SuperRunName $MultiFQFile 
	# [SaveDir] [ReferenceName] [SuperRunName] [ListOfFastqFiles]

fi


# Efterfølg med Annotation_report.rmd