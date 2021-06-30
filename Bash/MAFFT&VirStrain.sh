
# Multiple sequence alignment med MAFFT
mafft-linsi 'in.fasta' > 'out.mafft' 
# mafft-linsi --clustalout 'in.fasta' > 'out.clustal' # For clustal format

conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain

# Ændrer stier og navne så de passer
MainF=/home/pato/Skrivebord/HPV16_projekt
VirStrain_customdb=$MainF/References/Combined_mainlines_wRevised_wHPV-mTypes
MAFFTIn=$MainF/References/0Andre/Combined_mainlines_wRevised_wHPV-mTypes.mafft

# Byggger database
python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1

#TEST
#FQin=/home/pato/Skrivebord/HPV16_projekt/FASTQ/Pt_11_DNA.IonXpress_085.fastq
#Resultsout=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/test_16substrain
#python VirStrain.py -i $FQin -d $VirStrain_customdb -o $Resultsout 

conda deactivate





# Bygge mange databaser (eks 1 for hver HPV type)
conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain
MainF=/home/pato/Skrivebord/HPV16_projekt

for hpvtype in $MainF/References/0Andre/HPV_all_types/MSA/HPV*.mafft; do

	# Ændrer stier og navne så de passer
	HPVName=${hpvtype##*/}
	HPVn=${HPVName%.mafft}
	VirDBName=$(echo ${HPVn}_VirStrainDB)
	VirStrain_customdb=$MainF/References/VirStrainDBs_HPV_sub/$VirDBName
	echo $VirDBName
	MAFFTIn=$hpvtype

	# Byggger database
	python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1

done

conda deactivate