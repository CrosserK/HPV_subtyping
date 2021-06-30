
# Finde alle som er maintypes og navngive gfffiler og deres chr name herefter

FilePath=/home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types


for file in $FilePath/All_important_HPV_MainAndSub_revised/*.txt; do

	#####TEST
	#file=$FilePath/All_important_HPV_MainAndSub_revised/HPV16.txt


	cFile=$file
	cHPV=${cFile##*/}
	# Samler alle fasta filer til en fil til mafft
	ListeMedIDs=$(cut -f5 $file)
	echo $ListeMedIDs
	# Finder filer

done

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

done

rm address_list.txt
rm address_list_fix.txt
