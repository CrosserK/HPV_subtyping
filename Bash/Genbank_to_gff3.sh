
bp_genbank2gff 'in.gb' --viral --stdout > 'out.gff3'
# Make sure the name is correct. The script seems to remove the ".1" from "K02718.1" in the gff3 file