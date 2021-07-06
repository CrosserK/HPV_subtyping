
###################### EDIT ######################
MainF=/home/pato/Skrivebord/HPV16_projekt
VirStrainMaindb=$MainF/References/Combined_mainlines_wRevised_wHPV-mTypes
SuperRunName=Karoline_92fastq_6revisedRefs_1416_02072021
##################################################


FastqList=$(< $MainF/FASTQfiles_${SuperRunName}.txt)

for f in $FastqList; do
	rm $MainF/VirStrain_run/$SuperRunName/${f}/Revised_SubTypeCalls.txt
	touch $MainF/VirStrain_run/$SuperRunName/${f}/Revised_SubTypeCalls.txt
	#echo "HPV56_X74483_1_revised" >> $MainF/VirStrain_run/$SuperRunName/${f}/Revised_SubTypeCalls.txt
done
