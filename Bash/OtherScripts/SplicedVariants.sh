
# Finder Splejsningsvarianter

# UD fra gff3 fil, find E4 splice koordinater

# Hent nukleotid info fra de to dele

OriNukl1="ATGCTAGCA"
OriNukl2="CATGATCGATGCATC"

# Lav ny amino syre ramme fra spliced gener

SpliceNukl= OriNukl1 + OriNukl2
# SpliceNukl="ATGCTAGCACATGATCGATGCATC"

# eg 820-880 + 3350-3800
Aminosyrer=Oversæt(SpliceNukl)

# Oversæt positioner der matcher i splicet gen til ny ramme

For variant in varianterIndenfor820-880; do
	# Original position
	variant
	#Ny position
	variant=variant


elif [ in del2nukl ]; then
	OriNukl2

