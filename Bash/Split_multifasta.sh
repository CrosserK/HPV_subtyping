
#Target folders:
MainF='/home/pato/Skrivebord/HPV16_projekt'
RefFolder=$MainF/References/GFFfiles/All_550

cd $RefFolder
cat multi.fasta | awk '{
    if (substr($0, 1, 17)==">") {filename=(substr($1,1) ".fasta")}
    print $0 > filename
}'


# Split multi gff3 fil
MultiGFF=All_550_HPV_norevised.gff3
grep "sequence-region" $MultiGFF | awk '{print $2}' > identifiers
for i in `cat ./identifiers`; do 
	grep -e $i -e Species -A1 $MultiGFF > $i.gff3
done

# Rename to remove anything after space
for f in *.fasta; do
	mv "$f" "$(echo $(echo $f | cut -f1 -d " " ).fasta)"
done

cd 

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



# Gør noget på alle filer i mappe (eks med at lave alle sekvenser i fasta til uppercase:

for i in *
do
    if test -f "$i" 
    then
    	awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' $i > ${i%.fasta}.uppercase.fasta	
    	echo "Making uppercase to $i"
    fi
done

# Rename all files ex:
for i in pt._*.fastq; do
	mv $i test/$(echo pt_${i#pt._})
done