
# Combine multiple lines of same gene because of "join". Get start position of first line and stop position of last line.


file='/home/pato/Skrivebord/HPV16_projekt/Andet/HPV70_U21941_1_mod.gff3' 


# Find those that have duplicates for genes (positions are different)

for file in /home/pato/Skrivebord/HPV16_projekt/References/GFFfiles/*.gff3; do
	
	Duplicates=$(awk '{print $9}' $file | uniq -d )
	if [ ${#Duplicates} -gt 0 ]; then
	echo THIS FILE HAS DUPS: $file
	echo "$Duplicates"
	fi

done


cat $file | awk 

# 



