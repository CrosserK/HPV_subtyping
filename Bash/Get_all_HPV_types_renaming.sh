

# NOTE: Der er 1 sted som ikke udføres automatisk, der skal laves lidt manuelt! (Der står hvor det er)

# Get all sublineages of a HPV supertype
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
File=$FilePath/PaVE_HPV_all_main.txt

for number in 16 18 31 33 35 39 45 51 52 56 58 59 68 26 53 66 67 70 73 82 30 34 69 85 97 6 11; do
	GETHPVType=$(echo HPV${number})
	grep -P "${GETHPVType}\t" $File > $FilePath/All_important_HPV_Main/${GETHPVType}.txt
done



# Remove all space
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types

for number in 16 18 31 33 35 39 45 51 52 56 58 59 68 26 53 66 67 70 73 82 30 34 69 85 97 6 11; do
	GETHPVType=$(echo HPV${number})
	File=$FilePath/All_important_HPV_Main/${GETHPVType}.txt
	sed 's/ \+//g' $File > $FilePath/All_important_HPV_Main_nospace/${GETHPVType}.txt
done

	#sed -E 's/[[:space:]]/-/g' $FilePath/${GETHPVType}.txt > $FilePath/All_important_HPV_Main/${GETHPVType}.txt


# Rearrange tabs
FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types

for number in 16 18 31 33 35 39 45 51 52 56 58 59 68 26 53 66 67 70 73 82 30 34 69 85 97 6 11; do
	GETHPVType=$(echo HPV${number})
	File=$FilePath/All_important_HPV_Main_nospace/${GETHPVType}.txt
	awk -F'\t' -v OFS="\t" '{ print $2, $1, $3, $5, $4, $6 }' $File > $FilePath/All_important_HPV_Main_fixed/${GETHPVType}.txt
done




FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types

# Sætter mainline øverst
for file in $FilePath/All_important_HPV_Main_fixed/*.txt; do
	cFile=$file
	cHPV=${cFile##*/}
	#echo $cHPV
	cat $file $FilePath/All_important_HPV_Sub/${cHPV} > temp && mv temp $FilePath/All_important_HPV_MainAndSub/${cHPV}
done

# Finder hvor mange der har revised
for file in $FilePath/All_important_HPV_Main_fixed/*.txt; do
	cFile=$file
	cHPV=${cFile##*/}
	grep "Revised" $FilePath/All_important_HPV_MainAndSub_revised/${cHPV}
done

# Indsætter de funde referencer manuelt og fjerner dem de erstatter
# Manuelt

FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types

# Multiple sequence alignment og variant calling på hver supertype med undertyper
for file in $FilePath/All_important_HPV_MainAndSub_revised/*.txt; do

	#####TEST
	#file=$FilePath/All_important_HPV_MainAndSub_revised/HPV16.txt


	cFile=$file
	cHPV=${cFile##*/}
	# Samler alle fasta filer til en fil til mafft
	ListeMedIDs=$(cut -f5 $file)
	#echo $ListeMedIDs
	# Finder filer
	rm address_list.txt
	for ID in $ListeMedIDs; do
	echo $FilePath/MainAndSub_merge_revisedInserted/$ID >> address_list.txt
	done

	adresses=$(< address_list.txt)

	rm address_list_fix.txt

	for element in $adresses; do
	echo "$element"\*.fasta >> address_list_fix.txt # Tilføjer stjerne og suffix
	done

	# Samler fasta filer til 1 fasta
	address_list_fix=$(< address_list_fix.txt)

	cat $address_list_fix > $FilePath/Samlet_fasta/${cHPV%txt}fasta 

	# cat har allerede fundet korrekte filer gennem * 
	# Begynder at finde faktisk navn på mainstrain e.g om de ender på .1 .2 eller uden nogen af dem, så det kan benyttes til chrname

	MainStrain=($address_list_fix) # Gemmer mainstrain navn, virker da den kun tager første instans
	RefBaseName=${MainStrain##*/} 
	ChrName=${RefBaseName%.fasta}


	# Aligner med mafft
	mafft-linsi $FilePath/Samlet_fasta/${cHPV%txt}fasta  > $FilePath/MSA/${cHPV%txt}mafft
	# Kalder varianter med snp-sites
	snp-sites -v -o $FilePath/SNPs/${cHPV%txt}vcf $FilePath/MSA/${cHPV%txt}mafft
	
	MainF=/home/pato/Skrivebord/HPV16_projekt
	Rscriptfolder=$MainF/Scripts/R
	mafftFile=$FilePath/MSA/${cHPV%txt}mafft
	mafftName=${mafftFile##*/}
	# Korrigerer
	Rscript $Rscriptfolder/Correct_MSA_vcf_File_to_main_strain_bash.R $FilePath $mafftName ${mafftName%mafft}vcf $ChrName

	CorrVCF=$FilePath/SNPs_corr/corr_${cHPV%txt}vcf

	# Omsætter til bed fil
	grep -v "#" $CorrVCF | awk '{print $1, $2, $2 + 1}' >  $FilePath/SNPs_as_bed/${cHPV%txt}bed

done
rm address_list.txt
rm address_list_fix.txt




FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types
MainF=/home/pato/Skrivebord/HPV16_projekt
# Gør ovenstående med ny navngivning (V5)

# Multiple sequence alignment og variant calling på hver supertype med undertyper
for file in $FilePath/All_important_HPV_MainAndSub/*.txt; do

	#####TEST
	#file=$FilePath/All_important_HPV_MainAndSub_revised/HPV16.txt

	cFile=$file
	cHPV=${cFile##*/}
	# Samler alle fasta filer til en fil til mafft
	ListeMedIDs=$(cut -f5 $file)
	#echo $ListeMedIDs
	# Finder subtype IDs som er under HPV overtype
	
	rm ID_list.txt
	for ID in $ListeMedIDs; do
	echo $ID >> ID_list.txt
	done

	IDs=$(< ID_list.txt)

	touch Empty.fasta
	
	NewFasta=Empty.fasta
	for f in $IDs; do
	cat $NewFasta $FilePath/V5_ALL_main_and_sub_revisedINS_rdy/*${f}*.fasta >> tempfasta.fasta
	mv tempfasta.fasta tempfasta1.fasta
	NewFasta=tempfasta1.fasta
	done
	mv $NewFasta $FilePath/Samlet_fasta_ny/${cHPV}.fasta

	#$FilePath/V5_ALL_main_and_sub_revisedINS_rdy/${cHPV%txt}fasta 

	#TEST
	# cat har allerede fundet korrekte filer gennem * 
	# Begynder at finde faktisk navn på mainstrain e.g om de ender på .1 .2 eller uden nogen af dem, så det kan benyttes til chrname
	# Henter første ID, som er main strain:
	MainStrain=$(awk 'sub(/^>/, "")' $FilePath/Samlet_fasta_ny/${cHPV}.fasta)
	arr=($MainStrain)
	ChrName=${arr[0]}

	# Aligner med mafft
	#mafft-linsi $FilePath/Samlet_fasta_ny/${cHPV}.fasta  > $FilePath/MSA/${cHPV%txt}mafft
	# Kalder varianter med snp-sites
	#snp-sites -v -o $FilePath/SNPs/${cHPV%txt}vcf $FilePath/MSA/${cHPV%txt}mafft
	
	MainF=/home/pato/Skrivebord/HPV16_projekt
	Rscriptfolder=$MainF/Scripts/R
	mafftFile=$FilePath/MSA/${cHPV%txt}mafft
	mafftName=${mafftFile##*/}
	# Korrigerer
	Rscript $Rscriptfolder/Correct_MSA_vcf_File_to_main_strain_bash.R $FilePath $mafftName ${mafftName%mafft}vcf $ChrName

	CorrVCF=$FilePath/SNPs_corr/corr_${cHPV%txt}vcf

	# Omsætter til bed fil
	grep -v "#" $CorrVCF | awk '{print $1, $2, $2 + 1}' >  $FilePath/SNPs_as_bed/${cHPV%txt}bed

	echo Færdig med MSA: $ListeMedIDs

done
rm Empty.fasta
rm ID_list.txt