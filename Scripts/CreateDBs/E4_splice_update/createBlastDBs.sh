

location=/home/pato/Skrivebord/E4_splice_update/Organised_Fastas

for folder in $location/*; do
	
	foldername=$(basename $folder)
	mkdir -p ${folder}/blastdb
	cat ${folder}/*.fasta > ${folder}/blastdb/${foldername}_db.fasta

	makeblastdb -dbtype nucl -in ${folder}/blastdb/${foldername}_db.fasta -title $foldername

done



# Run query to db
location=/home/pato/Skrivebord/E4_splice_update/Organised_Fastas
foldername=HPV68
queryfile=${location}/${foldername}/HPV68_E4_seq.fasta
db=/home/pato/Skrivebord/E4_splice_update/Organised_Fastas/HPV-mwg1c09/blastdb/HPV68_db.fasta


${db}
blastn -query $queryfile -db HPV68_db -task blastn-short

