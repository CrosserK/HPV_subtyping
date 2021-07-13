
# Finde alle som er maintypes og navngive gfffiler og deres chr navn herefter

FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types


for file in $FilePath/All_important_HPV_MainAndSub_revised/*.txt; do

	#####TEST
	file=$FilePath/All_important_HPV_MainAndSub_revised/HPV16.txt

	cFile=$file
	cHPV=${cFile##*/}
	# Samler alle fasta filer til en fil til mafft
	ListeMedIDs=$(cut -f5 $file)
	#echo $ListeMedIDs
	# Finder subtype IDs som er under HPV overtype

	for f in ${FilePath}/*.gff3; do 
	
	FileName=${f##*/}
	ChrName=${FileName%.gff3}
	NewChrName=$( sed 's/\./_/g' <<< "$ChrName" ) # <<< is a here-string
	
	# Ændrer chrnavn i fil
	sed -i "s/$ChrName/$NewChrName/g" $f

	# Ændrer navn
	mv $f $FilePath/${NewChrName}.gff3
	# rename -n 's/\./_/g' $f

	done
	
	#Finder GenBank ID for maintyper

	# Henter første ID, som er main strain:
	MainStrain=$(awk 'sub(/^>/, "")' $FilePath/Samlet_fasta_ny/${cHPV}.fasta)
	arr=($MainStrain)
	ChrName=${arr[0]}

	NewName=$(cat $NewFasta $FilePath/V5_ALL_main_and_sub_revisedINS_rdy/*${ChrName}*.fasta | awk 'FNR == 1 {print $1}' | awk 'sub(/^>/, "")')

	sed -i "s/$ChrName/$NewChrName/g" $f

	mv $NewFasta $FilePath/Samlet_fasta_ny/${cHPV}.fasta

	#$FilePath/V5_ALL_main_and_sub_revisedINS_rdy/${cHPV%txt}fasta 

	#TEST
	# cat har allerede fundet korrekte filer gennem * 
	# Begynder at finde faktisk navn på mainstrain e.g om de ender på .1 .2 eller uden nogen af dem, så det kan benyttes til chrname
	# Henter første ID, som er main strain:
	MainStrain=$(awk 'sub(/^>/, "")' $FilePath/Samlet_fasta_ny/${cHPV}.fasta)
	arr=($MainStrain)
	ChrName=${arr[0]}

	echo Færdig med MSA: $ListeMedIDs

done

rm address_list.txt
rm address_list_fix.txt
