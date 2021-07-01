

# Finder HPV maintyper ud fra tabel fra PaVE. Hver maintype er tilknyttet et genbank ID (eg. K02718), men det fulde ID står ikke i tabellen. Der mangler ofte ".1" 
# eller ".2". Navngiver alle hovedtyper fra GenbankID til HPVx_GenbankID_y, hvor x er HPV nummer og y er nummeret for ".1" eller ".2" og gør det samme for
# deres chrname i filen. For dem som har revised PaVE versioner, tilføjes også "_revised" til endelsen. Dette gøres for både .fasta og .gff3 filer
# Har en tabel med HPV main typer fra PaVE og har en mappe med alle fasta'er og gff3'er navngivet med GenbankID 


# Get all mainlines of a HPV supertype from text table
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
File=$FilePath/Pave_MainLines_compact.txt



# Ændrer fasta "_" til "."
for f in $FilePath/PaVE_Main\&Sub/*.fasta; do 
	
	FileName=${f##*/}
	ChrName=${FileName%.fasta}
	NewChrName=$( sed 's/\./_/g' <<< "$ChrName" ) # <<< is a here-string
	
	# Ændrer chrnavn i fil
	sed -i "s/$ChrName/$NewChrName/g" $f

	# Ændrer navn
	mv $f $FilePath/PaVE_Main\&Sub/${NewChrName}.fasta
	# rename -n 's/\./_/g' $f

done

# Fjerner alt efter GenbankID i fasta header
for f in $FilePath/PaVE_Main\&Sub/*.fasta; do 
	
	sed -i '/^>/ s/ .*//' $f

done

# Ændrer gff3 "_" til "."
for f in $FilePath/PaVE_Main\&SubGFF3/*.gff3; do 
	
	FileName=${f##*/}
	ChrName=${FileName%.gff3}
	NewChrName=$( sed 's/\./_/g' <<< "$ChrName" ) # <<< is a here-string
	
	# Ændrer chrnavn i fil
	sed -i "s/$ChrName/$NewChrName/g" $f

	# Ændrer navn
	mv $f $FilePath/PaVE_Main\&SubGFF3/${NewChrName}.gff3

done



# Ændrer nu alle filer sådan at hvis de er en "main" type, kommer de på formen HPVx_GenbankID_y, både i navn og i fil


START=2
END=$(awk 'END{print NR}' $File)
for (( linenumber=$START; linenumber<=$END; linenumber++ )); do

	# Finder første HPV main nummer, tilhørende genbank ID og om PaVE har en revised version
	HPVName=$(cat $File | awk -v awkvar=$linenumber 'FNR == awkvar {print $1}')
	GenBankIDmiss=$(cat $File | awk -v awkvar=$linenumber 'FNR == awkvar {print $2}')
	Revised=$(cat $File | awk -v awkvar=$linenumber 'FNR == awkvar {print $3}')

	#echo $HPVName $GenBankIDmiss $Revised






	# Ændrer filnavn og chrname i fasta'er
	RefF=$FilePath/PaVE_Main\&Sub

	FFile=$(find $RefF -name ${GenBankIDmiss}*.fasta -type f)
	Fdir=$(dirname $FFile)
	FF=${FFile##*/} # Finder basename
	RealID=${FF%.fasta}
	
	# Laver "."" til "_"
	RealID_nodot=$(sed 's/\./_/g' <<< $RealID) # <<< gives a here-string

	if [ "$Revised" == "Revised" ]; then

	# Ændrer filnavn
	mv $FFile $Fdir/${HPVName}_${RealID_nodot}_revised.fasta

	# Ændrer chrnavn i fil til HPVx_GenbankID_y
	RealName=${HPVName}_${RealID_nodot}_revised

	sed -i "s/$RealID/$RealName/g" $Fdir/${HPVName}_${RealID_nodot}_revised.fasta
	echo IMPORTANT: Content of $HPVName fasta must be manually edited to contain revised version

	else

	# Ændrer filnavn
	mv $FFile $Fdir/${HPVName}_${RealID_nodot}.fasta

	# Ændrer chrnavn i fil til HPVx_genbankid_y
	RealName=${HPVName}_${RealID_nodot}

	sed -i "s/$RealID/$RealName/g" $Fdir/${HPVName}_${RealID_nodot}.fasta
 
	fi






	# Ændrer filnavn og chrname i gff3'er
	RefF=$FilePath/PaVE_Main\&SubGFF3

	FFile=$(find $RefF -name ${GenBankIDmiss}*.gff3 -type f)
	Fdir=$(dirname $FFile)
	FF=${FFile##*/} # Finder basename
	RealID=${FF%.gff3}
	
	# Laver "."" til "_"
	RealID_nodot=$(sed 's/\./_/g' <<< $RealID) # <<< gives a here-string

	if [ "$Revised" == "Revised" ]; then

	# Ændrer filnavn
	mv $FFile $Fdir/${HPVName}_${RealID_nodot}_revised.gff3

	# Ændrer chrnavn i fil til HPVx_GenbankID_y
	RealName=${HPVName}_${RealID_nodot}_revised

	sed -i "s/$RealID/$RealName/g" $Fdir/${HPVName}_${RealID_nodot}_revised.gff3
	echo IMPORTANT: Content of $HPVName gff3 must be manually edited to contain revised version

	else

	# Ændrer filnavn
	mv $FFile $Fdir/${HPVName}_${RealID_nodot}.gff3

	# Ændrer chrnavn i fil til HPVx_genbankid_y
	RealName=${HPVName}_${RealID_nodot}

	sed -i "s/$RealID/$RealName/g" $Fdir/${HPVName}_${RealID_nodot}.gff3
 
	fi


done

# Der nogle i listen som giver en fejl, men det er fordi PaVE ikke har oplyst noget genbank ID, og der derfor ikke er nogen fasta eller gff3 fil tilgængelig

# Scriptet outputtede de revised filer og de skal selvfølgelig opdateres. Det gøres manuelt ved at hente dem fra PaVE.




# Laver nu multiple sequence alignment på subtyper efter maintype
# Grupperer efter maintyper

SubLineFile=$FilePath/PaVE_HPV_all_sublineages.txt
MainLineFile=$FilePath/Pave_MainLines_compact.txt
RefF=$FilePath/PaVE_Main\&Sub

mkdir -p $FilePath/Combined_Fastas
mkdir -p $FilePath/MSAs

START=2
END=$(awk 'END{print NR}' $MainLineFile)
for (( linenumber=$START; linenumber<=$END; linenumber++ )); do


	# Finder første HPV main nummer, tilhørende genbank ID og om PaVE har en revised version
	HPVName=$(cat $MainLineFile | awk -v awkvar=$linenumber 'FNR == awkvar {print $1}')
	GenBankIDmiss=$(cat $MainLineFile | awk -v awkvar=$linenumber 'FNR == awkvar {print $2}')
	Revised=$(cat $MainLineFile | awk -v awkvar=$linenumber 'FNR == awkvar {print $3}')

	#TEST
	#HPVName=HPV1

	echo $HPVName

	# Finder alle subtyper under maintype
	AllTypesList=($(cat $SubLineFile | grep "$HPVName\s" | awk '{print $5}')) # \s er for space. Sikrer at eks HPV16 ikke matcher HPV160


	declare -a AllAddressList=($(find $RefF -name "${HPVName}_*" -type f))

	if [ ${#AllAddressList} -gt 0 ]; then # Tester om der er fundet en fil (Der var en genbank ID ved det HPV nummer)
	# Finder de pågældende filer

	# Laver liste over adresser som skal gives til cat.  " " er for at give mellemrum mellem appending til liste
	for f in "${AllTypesList[@]}"; do
	AllAddressList+=" "
	AllAddressList+=$( find $RefF -name "*${f}*" -type f )
	done

	cat $AllAddressList > $FilePath/Combined_Fastas/${HPVName}.fasta

	# Fjerner evt. duplikater der kan være opstået fordi nogen main refs også er i subtyper refs
	awk '/^>/{f=!d[$1];d[$1]=1}f' $FilePath/Combined_Fastas/${HPVName}.fasta > $FilePath/Combined_Fastas/${HPVName}.tmp
	mv $FilePath/Combined_Fastas/${HPVName}.tmp $FilePath/Combined_Fastas/${HPVName}.fasta

	# Fjerner blanke linjer
	sed -i '/^$/d' $FilePath/Combined_Fastas/${HPVName}.fasta

	# Fjern alle RYSWKMBDHV fra nukleotid sekvenser
	# sed -e '/^>/! s/[RYSWKMBDHV]/N/g' in.fasta > out.fasta 

	# Tester om der er mere en 1 reference i samlet fil. Ellers er der ingen undertyper og den skal ikke behandles af mafft

	MoreThanOne=$(grep ">" $FilePath/Combined_Fastas/${HPVName}.fasta | grep "" -c)

	if [ $MoreThanOne -gt 1 ]; then

	# Laver multiple sequence alignment med mafft
	mafft-linsi $FilePath/Combined_Fastas/${HPVName}.fasta > $FilePath/MSAs/${HPVName}.mafft
	
	# Kalder varianter med snp-sites
	#snp-sites -v -o $FilePath/SNPs/${cHPV%txt}vcf $FilePath/MSA/${cHPV%txt}mafft
	else	
	echo No subtypes for $HPVName

	fi


	fi


	# Laver VirStrain database fra subtyperede filer





done


# Laver VirStrain database for alle Main typer

find ^HPV





















#### GAMMELT

## NOTE: Der er 1 sted der skal udføres manuelt. Der står hvor det er.
#
## Get all sublineages of a HPV supertype from text table
#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
#File=$FilePath/PaVE_HPV_all_main.txt
#HPVlist="16 18 31 33 35 39 45 51 52 56 58 59 68 26 53 66 67 70 73 82 30 34 69 85 97 6 11"
#
#
#for number in $HPVlist; do
#	GETHPVType=$(echo HPV${number})
#	grep -P "${GETHPVType}\t" $File > $FilePath/All_important_HPV_Main/${GETHPVType}.txt
#done
#
#
#
#
## Remove all spaces
#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
#
#for number in $HPVlist; do
#	GETHPVType=$(echo HPV${number})
#	File=$FilePath/All_important_HPV_Main/${GETHPVType}.txt
#	sed 's/ \+//g' $File > $FilePath/All_important_HPV_Main_nospace/${GETHPVType}.txt
#done
#
#
#
#
#
## Rearrange tabs
#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
#
#for number in $HPVlist; do
#	GETHPVType=$(echo HPV${number})
#	File=$FilePath/All_important_HPV_Main_nospace/${GETHPVType}.txt
#	awk -F'\t' -v OFS="\t" '{ print $2, $1, $3, $5, $4, $6 }' $File > $FilePath/All_important_HPV_Main_fixed/${GETHPVType}.txt
#done
#
#
#
#
#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
#
## Sætter mainline øverst
#for file in $FilePath/All_important_HPV_Main_fixed/*.txt; do
#	cFile=$file
#	cHPV=${cFile##*/}
#	#echo $cHPV
#	cat $file $FilePath/All_important_HPV_Sub/${cHPV} > temp && mv temp $FilePath/All_important_HPV_MainAndSub/${cHPV}
#done
#
#
#
#
## Finder hvor mange der har revised
#for file in $FilePath/All_important_HPV_Main_fixed/*.txt; do
#	cFile=$file
#	cHPV=${cFile##*/}
#	grep "Revised" $FilePath/All_important_HPV_MainAndSub_revised/${cHPV}
#done
#
#
#
#
#
#
## Indsæt de funde referencer manuelt og fjern dem de erstatter
#
#
#
## Multiple sequence alignment og variant calling på hver supertype med undertyper
#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
#MainF=/home/pato/Skrivebord/HPV16_projekt
## Gør ovenstående med ny navngivning (V5)
#
#for file in $FilePath/All_important_HPV_MainAndSub/*.txt; do
#
#	#####TEST
#	#file=$FilePath/All_important_HPV_MainAndSub_revised/HPV16.txt
#
#	cFile=$file
#	cHPV=${cFile##*/}
#	# Samler alle fasta filer til en fil til mafft
#	ListeMedIDs=$(cut -f5 $file)
#	#echo $ListeMedIDs
#	# Finder subtype IDs som er under HPV overtype
#	
#	rm ID_list.txt
#	for ID in $ListeMedIDs; do
#	echo $ID >> ID_list.txt
#	done
#
#	IDs=$(< ID_list.txt)
#
#	touch Empty.fasta
#	
#	NewFasta=Empty.fasta
#	#Finder fasta fra GenBankID
#	for f in $IDs; do
#	cat $NewFasta $FilePath/V5_ALL_main_and_sub_revisedINS_rdy/*${f}*.fasta >> tempfasta.fasta
#	mv tempfasta.fasta tempfasta1.fasta
#	NewFasta=tempfasta1.fasta
#	done
#	mv $NewFasta $FilePath/Samlet_fasta_ny/${cHPV}.fasta
#
#	#$FilePath/V5_ALL_main_and_sub_revisedINS_rdy/${cHPV%txt}fasta 
#
#	#TEST
#	# cat har allerede fundet korrekte filer gennem * 
#	# Begynder at finde faktisk navn på mainstrain e.g om de ender på .1 .2 eller uden nogen af dem, så det kan benyttes til chrname
#	# Henter første ID, som er main strain:
#	MainStrain=$(awk 'sub(/^>/, "")' $FilePath/Samlet_fasta_ny/${cHPV}.fasta)
#	arr=($MainStrain)
#	ChrName=${arr[0]}
#
#	# Aligner med mafft
#	#mafft-linsi $FilePath/Samlet_fasta_ny/${cHPV}.fasta  > $FilePath/MSA/${cHPV%txt}mafft
#	# Kalder varianter med snp-sites
#	#snp-sites -v -o $FilePath/SNPs/${cHPV%txt}vcf $FilePath/MSA/${cHPV%txt}mafft
#	
#	MainF=/home/pato/Skrivebord/HPV16_projekt
#	Rscriptfolder=$MainF/Scripts/R
#	mafftFile=$FilePath/MSA/${cHPV%txt}mafft
#	mafftName=${mafftFile##*/}
#	# Korrigerer
#	Rscript $Rscriptfolder/Correct_MSA_vcf_File_to_main_strain_bash.R $FilePath $mafftName ${mafftName%mafft}vcf $ChrName
#
#	CorrVCF=$FilePath/SNPs_corr/corr_${cHPV%txt}vcf
#
#	# Omsætter til bed fil
#	grep -v "#" $CorrVCF | awk '{print $1, $2, $2 + 1}' >  $FilePath/SNPs_as_bed/${cHPV%txt}bed
#
#	echo Færdig med MSA: $ListeMedIDs
#
#done
#
## Fjerner duplicate referencer, de kan være opstået fordi der er nogen maintyper i tabel med undertypevarianter
#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/MSA
#
#for f in $FilePath/*.mafft; do
#	awk '/^>/{f=!d[$1];d[$1]=1}f' $f > ${f}.tmp
#	mv ${f}.tmp $f
#done
#
#
#
#rm Empty.fasta
#rm ID_list.txt















# GAMMELT

#FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
#
## Multiple sequence alignment og variant calling på hver supertype med undertyper
#for file in $FilePath/All_important_HPV_MainAndSub_revised/*.txt; do
#
#	#####TEST
#	#file=$FilePath/All_important_HPV_MainAndSub_revised/HPV16.txt
#
#
#	cFile=$file
#	cHPV=${cFile##*/}
#	# Samler alle fasta filer til en fil til mafft
#	ListeMedIDs=$(cut -f5 $file)
#	#echo $ListeMedIDs
#	# Finder filer
#	rm address_list.txt
#	for ID in $ListeMedIDs; do
#	echo $FilePath/MainAndSub_merge_revisedInserted/$ID >> address_list.txt
#	done
#
#	adresses=$(< address_list.txt)
#
#	rm address_list_fix.txt
#
#	for element in $adresses; do
#	echo "$element"\*.fasta >> address_list_fix.txt # Tilføjer stjerne og suffix
#	done
#
#	# Samler fasta filer til 1 fasta
#	address_list_fix=$(< address_list_fix.txt)
#
#	cat $address_list_fix > $FilePath/Samlet_fasta/${cHPV%txt}fasta 
#
#	# cat har allerede fundet korrekte filer gennem * 
#	# Begynder at finde faktisk navn på mainstrain e.g om de ender på .1 .2 eller uden nogen af dem, så det kan benyttes til chrname
#
#	MainStrain=($address_list_fix) # Gemmer mainstrain navn, virker da den kun tager første instans
#	RefBaseName=${MainStrain##*/} 
#	ChrName=${RefBaseName%.fasta}
#
#
#	# Aligner med mafft
#	mafft-linsi $FilePath/Samlet_fasta/${cHPV%txt}fasta  > $FilePath/MSA/${cHPV%txt}mafft
#	# Kalder varianter med snp-sites
#	snp-sites -v -o $FilePath/SNPs/${cHPV%txt}vcf $FilePath/MSA/${cHPV%txt}mafft
#	
#	MainF=/home/pato/Skrivebord/HPV16_projekt
#	Rscriptfolder=$MainF/Scripts/R
#	mafftFile=$FilePath/MSA/${cHPV%txt}mafft
#	mafftName=${mafftFile##*/}
#	# Korrigerer
#	Rscript $Rscriptfolder/Correct_MSA_vcf_File_to_main_strain_bash.R $FilePath $mafftName ${mafftName%mafft}vcf $ChrName
#
#	CorrVCF=$FilePath/SNPs_corr/corr_${cHPV%txt}vcf
#
#	# Omsætter til bed fil
#	grep -v "#" $CorrVCF | awk '{print $1, $2, $2 + 1}' >  $FilePath/SNPs_as_bed/${cHPV%txt}bed
#
#done
#rm address_list.txt
#rm address_list_fix.txt
