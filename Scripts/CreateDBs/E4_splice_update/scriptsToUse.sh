
# For each subtype
for gff in *.gff3; do
	if [[ "$gff" != *"HPV"* ]];then
    
	subtype=${gff%.gff3}

	# Find maintype name
	echo looking for $subtype

	path=$(find /home/pato/Skrivebord/HPV_subtyping/Scripts/CreateDBs/E4_splice_update/Organised_Fastas/ -maxdepth 2 -name ${subtype}.fasta)
	dir=$(dirname $path)
	maintype=$(basename $dir)

	echo $maintype $subtype
	python /home/pato/Skrivebord/HPV_subtyping/Scripts/CreateDBs/E4_splice_update/createGffForSubtype.py -s $subtype -m $maintype

  fi
 done

# Delete inserted lines:
for file in *.gff3; do
	awk '!/AUHpred/' $file > tmpfile.csv && mv tmpfile.csv $file
done

# Remove previous E4 lines
for file in *.gff3; do
	if [[ "$gff" != *"HPV"* ]];then
	awk '!/-E4-/' $file > tmpfile.csv && mv tmpfile.csv $file
	fi
done
