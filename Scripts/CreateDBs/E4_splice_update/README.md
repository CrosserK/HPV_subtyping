
# Preperation
Use seqGet to get seqeuences from PaVE
Use gff3Get.ipynb to get gff3 from PaVE
Organize into subfolders with Organize_fastas

# Find E4 seq for subtypes
Use createE4seqFromMain_andAlign.py to get sequences of the mainline and align them to all subtypes. The coordinates must be entered manually in the script
Use createE4seqFromMain to fetch the E4 coordinates from the maintype.
Use to align against the subtype and create coordinates

# Update gffs
createGFFForSubtype will use the created E4FixCoords.csv to update the original gff3 and create the E4fix gff3 accordingly. 
s

