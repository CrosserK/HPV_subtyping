

# Finder HPV maintyper ud fra tabel fra PaVE. Hver maintype er tilknyttet et genbank ID (eg. K02718), men det fulde ID står ikke i tabellen. Der mangler ofte ".1" 
# eller ".2". Navngiver alle hovedtyper fra GenbankID til HPVx_GenbankID_y, hvor x er HPV nummer og y er nummeret for ".1" eller ".2" og gør det samme for
# deres chrname i filen. For dem som har revised PaVE versioner, tilføjes også "_revised" til endelsen. Dette gøres for både .fasta og .gff3 filer
# Har en tabel med HPV main typer fra PaVE og har en mappe med alle fasta'er i $FilePath/PaVE_Main&Sub og gff3'er i mappedn $FilePath/PaVE_Main&SubGFF3 
# begge filformater navngivet efter GenbankID (e.g. "K02718.1.fasta")


# Get all mainlines of a HPV supertype from text table
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
File=$FilePath/Pave_MainLines_compact.txt



# Ændrer fasta "." til "_"
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


# Laver VirStrain database for alle subtyper. 
conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain
mkdir -p $FilePath/SubTypesCombined

# Gør det først for alle HPV som har et nummer og derefter for alle HPV som har formatet HPV-mxxxxx
for f in $FilePath/Combined_Fastas/HPV[0-9]*.fasta; do
	
	# Finder basename
	HPVName=${f##*/}
	mafft-linsi $f > $FilePath/MSAs/${HPVName%fasta}mafft
	# Byggger database
	VirStrain_customdb=$MainF/References/VirStrainDBs_HPV_sub/${HPVName%.fasta}_VirStrainDB
	MAFFTIn=$FilePath/MSAs/${HPVName%fasta}mafft
	python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1

done

for f in $FilePath/Combined_Fastas/HPV-m*.fasta; do
	
	# Finder basename
	HPVName=${f##*/}
	mafft-linsi $f > $FilePath/MSAs/${HPVName%fasta}mafft
	# Byggger database
	VirStrain_customdb=$MainF/References/VirStrainDBs_HPV_sub/${HPVName%.fasta}_VirStrainDB
	MAFFTIn=$FilePath/MSAs/${HPVName%fasta}mafft
	python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1

done


# Laver VirStrain database for alle Main typer
RefF=$FilePath/PaVE_Main\&Sub
mkdir -p $FilePath/MainsCombined
# Samler 
cat $RefF/HPV*.fasta > $FilePath/MainsCombined/AllMainCombined_wRevised_wHPVm-types.fasta
# MSA
mafft-linsi $FilePath/MainsCombined/AllMainCombined_wRevised_wHPVm-types.fasta > $FilePath/MainsCombined/AllMainCombined_wRevised_wHPVm-types.mafft
# Ændrer stier og navne så de passer
MainF=/home/pato/Skrivebord/HPV16_projekt
VirStrain_customdb=$MainF/References/AllMainCombined_wRevised_wHPVm-types
MAFFTIn=$FilePath/MainsCombined/AllMainCombined_wRevised_wHPVm-types.mafft
# Byggger database
python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1
conda deactivate





# Fikser revised GFF3 filer. De revised filer findes kun som genbank fil fra PaVE. 
# De skal derfor hentes ned manuelt og ligges i $FilePath/Revised_Genbank_files, med korrekte navne som script ovenover har givet, men med endelsen .gb.
# De kan da omdannes fra genbank til gff3 filer med følgende script. 
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types

for f in $FilePath/Revised_Genbank_files/*.gb; do 

	BaseName=${f##*/}
	BaseN=${BaseName%.gb}

	bp_genbank2gff $f --viral --stdout > ${f%gb}gff3
	# The script seems to remove the ".1" from "K02718.1" in the gff3 file

	# Læser HPVnummer
	HPVName=$(echo "$BaseN" | grep -o "HPV[0-9]*")

	# HUSK følgende 3 korrigeringer!
	# Lav korrekt navn, eks:
	sed -i "s/${HPVName}REF/$BaseN/g" ${f%gb}gff3

	# Fjern .t00:
	sed -i 's/.t00//g' ${f%gb}gff3

	# Fjern "_" fra gennavne for at annotation R script ikke kommer til at separere forkert. Der må gerne være "-"" i gennavne, men ikke i referencenavne

	# Sætert CDS phases til 0. De må ikke være "."
	# 1 tallet gør at den redigerer i det valgte, men også returnerer resten som ikke er redigeret
	awk '$3=="CDS" {print $8="0"} 1' OFS="\t" ${f%gb}gff3 > tmp.gff && mv tmp.gff ${f%gb}gff3
	# Fjerner enkelte linjer der opstår hvor der kun står "0"
	grep -v "^0" ${f%gb}gff3 > tmp.gff && mv tmp.gff ${f%gb}gff3
	# Fjerner nukleotid sekvens + dens header
	cat ${f%gb}gff3 | sed -n '/[>]/q;p' > tmp.txt 	# This says "when you reach the line that matches the pattern quit, otherwise print each line". The -n option prevents implicit printing and the p command is required to explicitly print lines.

	# Flytter klar fil til anden mappe
	mv tmp.txt $FilePath/Revised_GFF3_files_rdy/${BaseN}.gff3


done


# Finder hvilke gff3 filer fra NCBI der mangler gene rækker (nogle har kun rækker for CDS)
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/GFFfiles
for f in $FilePath/*.gff3; do 

	MustBeFixed=$(awk '{print $3}' $f | grep "gene")

	BaseName=${f##*/}

	if [ ${#MustBeFixed} -lt 1 ]; then

	echo $BaseName must be updated with genes

	fi

done


# Kopierer manuelt gb filer ned fra PaVE og benytter til følgende


# For genbank filer der skal laves om til gff3 filer
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/GFFfiles
for f in $FilePath/gbFilesToUpdateGff3/*.gb; do

	bp_genbank2gff $f --viral --stdout > ${f%gb}gff3 

	BaseName=${f##*/}
	BaseName=${BaseName%.gb}

	# Ændrer chr navn i fil fra NCBIs standard HPV[0-9]*REF, eg: HPV16REF til format som filerne er navngivet i
	sed -i "s/HPV[0-9]*REF/$BaseName/g" ${f%gb}gff3

	# Fjern .t00:
	sed -i 's/.t00//g' ${f%gb}gff3

	# Sætert CDS phases til 0. De må ikke være "."
	# 1 tallet gør at den redigerer i det valgte, men også returnerer resten som ikke er redigeret
	awk '$3=="CDS" {print $8="0"} 1' OFS="\t" ${f%gb}gff3 > tmp.gff && mv tmp.gff ${f%gb}gff3
	# Fjerner enkelte linjer der opstår hvor der kun står "0"
	grep -v "^0" ${f%gb}gff3 > tmp.gff && mv tmp.gff ${f%gb}gff3
	# Fjerner nukleotid sekvens + dens header
	cat ${f%gb}gff3 | sed -n '/[>]/q;p' > tmp.txt 	# This says "when you reach the line that matches the pattern quit, otherwise print each line". The -n option prevents implicit printing and the p command is required to explicitly print lines.
	mv tmp.txt ${f%gb}gff3

done