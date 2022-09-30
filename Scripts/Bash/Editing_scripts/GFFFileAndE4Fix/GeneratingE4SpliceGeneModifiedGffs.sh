
# Script to do most of the work in creating modified gff3 E4 splice gen fixes
# Put original gff in gffs folder before start

MainF=/home/pato/Skrivebord/HPV16_projekt
Folder=$MainF/gffs


for f in $Folder/*.gff3;do 
	
	basen=$(basename $f)
	hpvname=${basen%.gff3}

	# Fjerner alle linjer der har E4 i sig
	sed -i '/E4/d' $f

	# Indsætter E4 splice template linjer
	cat $f $MainF/GFFMAPPEtemplate.gff3 > ${f}.mod
	sed -i "s/HPV16_K02718_1_revised/$hpvname/g" ${f}.mod
	mv ${f}.mod $f

	# Laver E4fix fil
	cat $MainF/E4FIXCOORDSMAPPE.gff3 > $MainF/References/E4fixCoords/$basen

	sed -i "s/HPV16_K02718_1_revised/$hpvname/g" $MainF/References/E4fixCoords/$basen	

done