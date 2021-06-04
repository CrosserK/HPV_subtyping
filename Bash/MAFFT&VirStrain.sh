
# Multiple sequence alignment med MAFFT
mafft-linsi 'in.fasta' > 'out.afa' 
mafft-linsi --clustalout 'in.fasta' > 'out.clustal' # For clustal format

conda activate VarStrain # VirStrain har meget specifikke krav til pakke versioner, derfor køres conda env
cd VirStrain
# python VirStrain_build.py -i '/home/pato/Skrivebord/HPV16_projekt/References_andre/15_HPV_sublineages.mafft' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_15_virstrain'
python VirStrain.py  -m -i '/home/pato/Skrivebord/HPV16_projekt/FASTQ/Pt_33_RNA.IonXpress_087.fastq' -d '/home/pato/Skrivebord/HPV16_projekt/References/HPV16_15_virstrain' -o '/home/pato/Skrivebord/HPV16_projekt/VirStrain_run/test_m_param1'
conda deactivate