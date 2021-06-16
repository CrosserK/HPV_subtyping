
#Target folders:
MainF='/home/pato/Skrivebord/HPV16_projekt'
RefFolder=$MainF/References/HPV_revised_PaVE/HPV_HG19_V4

cd $RefFolder
cat HPV_HG19_V4_1.fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 > filename
}'


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
