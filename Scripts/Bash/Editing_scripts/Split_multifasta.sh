
#Target folders:
MainF='/home/pato/Skrivebord/HPV16_projekt'
RefFolder=$MainF/References/0Andre/HPV_all_types/PaVE_HPV_all_mainlines_split

MainF=/home/pato/Skrivebord/E4_splice_update/Combined_Fastas
RefFolder=/home/pato/Skrivebord/E4_splice_update/Organised_Fastas

# put all multifastas into a directory each
for x in ./*.fasta; do
  mkdir "${x%.*}" && mv "$x" "${x%.*}"
done

# Then split all to individual fastas
for HPV in *; do
cd $HPV
cat ${HPV}.fasta | awk '{ 
if (substr($0, 1, 1)==">") {filename=(substr($1,1) ".fasta")}
print $0 > filename
}'
rm ${HPV}.fasta
for i in *.fasta; do
	mv $i ${i#>} # Fjerner > fra navn
	mv $i ${HPV}_${i#>}
done
cd ..
done

mkdir HPV1_all




# Split multi gff3 fil
MultiGFF=All_550_HPV_norevised.gff3
grep "sequence-region" $MultiGFF | awk '{print $2}' > identifiers
for i in `cat ./identifiers`; do 
	grep -e $i -e Species -A1 $MultiGFF > $i.gff3
done

cd /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/PaVE_HPV_all_sublineages_split_renamed

# Rename to remove anything after space and do same for header in file
for f in *.fasta; do
	NewName="$(echo $(echo $f | cut -f1 -d " " ))"
	HeadName=$( echo \>${NewName%.fasta})
	# rename header
	echo $HeadName $NewName
	sed -i "1s/.*/$HeadName/" $f
done

cd /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/PaVE_HPV_all_sublineages_split_renamed

# Change name to remove anything before and after "_" including "_" and remove those files in another folder:
# Remove files containing:
for f in *.fasta; do
	# Gør det i to omgange fordi nogen har flere "_" 
	Startrm="${f%_*}"
	Endrm="${Startrm#*_}"
	Startrm2="${Endrm%_*}"
	Endrm2="${Endrm#*_}"
	# Gør det en gang for .1 og en gang for .2 endelser
	rm -f ../Main_\&_sub_merge_revisedInserted/"$(echo "$Endrm2".1.fasta)"
	rm -f ../Main_\&_sub_merge_revisedInserted/"$(echo "$Endrm2".2.fasta)"
done




# Combine fasta:

cat *fasta > combined.fasta



# Gør noget på alle filer i mappe (eks med at lave alle sekvenser i fasta til uppercase):

for i in *
do
    if test -f "$i" 
    then
    	awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' $i > ${i%.fasta}.uppercase.fasta
    	mv ${i%.fasta}.uppercase.fasta $i
    	echo "Making $i uppercase"
    fi
done

# Rename all files ex:
for i in pt._*.fastq; do
	mv $i test/$(echo pt_${i#pt._})
done

# Fjerner > fra navn
for i in *.fasta; do
	mv $i ${i#>} # Fjerner > fra navn
done

# Rename så . bliver til _ og ændrer header på samme måde
for f in *.fasta; do 
	Name=${f##*/}
	AName=${Name%.fasta}
	cat $f | sed 's/\./_/g' > ${f}.tmp
	mv ${f}.tmp $f 

	#Renaming
	i="${f%.fasta}"
	mv -i -- "$f" "${i//./_}.fasta"

done




# Rename med HPV nummer fra tabel
for f in *.fasta; do
	RefName=${f##*/}
	Shortname=${RefName%_*}
	# echo $Shortname

	SuperType=$(grep "$Shortname" ../download_human_RefClone_95d89352.txt | cut -f1)

	mv $f ${SuperType}_$RefName
done

cd /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/Revised_split

# Remame revised files
for f in *.fasta; do
	HPVType=${f%REF*}
	Startrm="${f%_*}"
	Endrm="${Startrm#*_}"
	Startrm2="${Endrm%_*}"
	Endrm2="${Endrm#*_}"
	NewName=$( echo ${HPVType}_${Endrm2}_1_revised.fasta)

	# Ændrer header
	HeadName=$( echo \>${NewName%.fasta})
	cp $f /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/Revised_split_renamed/$NewName
	sed -i "1s/.*/$HeadName/" /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/Revised_split_renamed/$NewName

done

# Derefter for special case:
for f in *-*.fasta; do
	HPVType=${f%KP692117*}
	NewName=$( echo ${HPVType}_${Endrm2}_1_revised.fasta)
	Startrm="${f%_*}"
	Endrm="${Startrm#*_}"
	Startrm2="${Endrm%_*}"
	Endrm2="${Endrm#*_}"
	NewName=$( echo ${HPVType}${Endrm2}_1_revised.fasta)

	# Ændrer header
	HeadName=$( echo \>${NewName%.fasta})
	cp $f /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/Revised_split_renamed/$NewName
	sed -i "1s/.*/$HeadName/" /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/Revised_split_renamed/$NewName

done

cd /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/PaVE_HPV_all_mainlines_split

# Fikser header i mainlines
for f in *.fasta; do
	RefName=${f##*/}
	RefN=$( echo \>${RefName%.fasta})
	sed -i "1s/.*/$RefN/" $f
done

# Fjerner dem som har en revised version fra submappe:
cd /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/V5_ALL_main_revisedINS
for f in *.fasta; do 
	RefName=${f##*/}
	NavnUdenRev=${RefName%_revised.fasta}
	# Finder GB ID
	Startrm="${NavnUdenRev%_*}"
	Endrm="${Startrm#*_}"
	Startrm2="${Endrm%_*}"
	Endrm2="${Endrm#*_}"
	rm ../V5_PaVE_HPV_all_sublineages_split_renamed/*${Endrm2}*.fasta
done

cd /home/pato/Skrivebord/HPV16_projekt/References/0Andre/HPV_all_types/PaVE_HPV_all_sublineages_split_renamed

# Flytter alle fra sublineage mappe til main_and_sub mappe
for f in *.fasta; do
	cp ../V5_PaVE_HPV_all_sublineages_split_renamed/$f ../V5_ALL_main_and_sub_revisedINS_rdy/$f
done

