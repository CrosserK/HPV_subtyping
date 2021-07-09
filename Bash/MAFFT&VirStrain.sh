
# Multiple sequence alignment med MAFFT
# For at lave MSA (multiple sequence alignment) bruges algoritmen mafft-linsi hvorefter man giver den 1 fasta fil med alle referencer der skal alignes til hinanden
# og peger med > hvor den alignede fil skal ligges og hvad den skal hedde. NOTE: ">" bruges til at henvise hvor resultater skal ligges, men man kan hurtig komme galt afsted
# og komme til at overskrive andre filer hvis man ikke passer på med hvad man kalder resultatfilen. Sørg altid for at navnet er noget som der ikke ligge i den sti man
# henviser til, da det ellers overskrives uden varsel. Eks på kørsel: 
mafft-linsi /path/to/in.fasta > another/path/to/out.mafft 
# mafft-linsi --clustalout /path/to/in.fasta > another/path/to/out.clustal # For clustal format

# Når den alignede fil er lavet, henvises VirStrain til denne fil ved at ændre i nedenstående adresse efter "MAFFTIn=".
# Herefeter ændres adressen og navnet efter "VirStrain_customdb=" til et passende sted og navn som kan gives til VirStrain i en subtypering
# Herefter kopieres alle linjer herunder indtil det første "conda deactivate" ind i terminalen (ctrl+shift+v, for at kopiere ind i terminal) og databasen bygges.

conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain

# Ændrer stier og navne så de passer
MainF=/home/pato/Skrivebord/HPV16_projekt
MAFFTIn=$MainF/References/0Andre/Combined_mainlines_wRevised_wHPV-mTypes.mafft
VirStrain_customdb=$MainF/References/Combined_mainlines_wRevised_wHPV-mTypes


# Byggger database
python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1

conda deactivate







#TEST
#FQin=/home/pato/Skrivebord/HPV16_projekt/FASTQ/Pt_11_DNA.IonXpress_085.fastq
#Resultsout=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/test_16substrain
#python VirStrain.py -i $FQin -d $VirStrain_customdb -o $Resultsout 


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