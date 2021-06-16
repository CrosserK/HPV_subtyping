
# Multiple sequence alignment med MAFFT
mafft-linsi 'in.fasta' > 'out.afa' 
# mafft-linsi --clustalout 'in.fasta' > 'out.clustal' # For clustal format

conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain
FQin=/home/pato/Skrivebord/HPV16_projekt/FASTQ/Pt_11_DNA.IonXpress_085.fastq
VirStrain_customdb=/home/pato/Skrivebord/HPV16_projekt/References/HPV16_16_virstrain
Resultsout=/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/test_16substrain

python VirStrain_build.py -i '/home/pato/Skrivebord/HPV16_projekt/References_andre/16substrain_HPV16_sublineages.mafft' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_16_virstrain' -s 1
python VirStrain.py  -m -i $FQin -d $VirStrain_customdb -o $Resultsout 

conda deactivate