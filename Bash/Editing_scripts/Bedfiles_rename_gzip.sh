
for f in /home/pato/Skrivebord/HPV16_projekt/References/BedFiles/Original/*.bed; do
	gffname=${f##*/}
	in=$f
	out=/home/pato/Skrivebord/HPV16_projekt/References/BedFiles/$gffname
	sed 's/K02718.1/K02718.1_revised/g' $in > $out 
done

for f in /home/pato/Skrivebord/HPV16_projekt/References/BedFiles/*.bed; do
	gzip $f
done