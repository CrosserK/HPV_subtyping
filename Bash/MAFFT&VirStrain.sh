
# Multiple sequence alignment med MAFFT
mafft-linsi 'in.fasta' > 'out.afa' 
# mafft-linsi --clustalout 'in.fasta' > 'out.clustal' # For clustal format

conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain
MainF=/home/pato/Skrivebord/HPV16_projekt
VirStrain_customdb=$MainF/References/HPV16_16_virstrain_revised
MAFFTIn=$MainF/References/HPV16_16subtypes_revised.mafft

python VirStrain_build.py -i $MAFFTIn -d $VirStrain_customdb -s 1


#FQin=/home/pato/Skrivebord/HPV16_projekt/FASTQ/Pt_11_DNA.IonXpress_085.fastq
#Resultsout=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/test_16substrain
#python VirStrain.py -i $FQin -d $VirStrain_customdb -o $Resultsout 

conda deactivate