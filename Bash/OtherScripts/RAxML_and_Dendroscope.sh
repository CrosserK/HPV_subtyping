



# For MSA fil til RAxMML:
mafft-linsi --thread 12 .fa > .mafft 


# RAxMML and dendroscope


MainF=/home/pato/Skrivebord/HPV16_projekt
RAxDir=$MainF/RAxML
DNA=$MainF/References/0Andre/Combined_mainlines_wRevised_wHPV-mTypes.mafft
Threads=T12
NumberOfTrees=100
OutDir=$RAxDir/HPVALLtest
mkdir -p $OutDir


raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# $NumberOfTrees -s $DNA -n $Threads -w $OutDir
mv $OutDir/RAxML_bipartitions.T12 $OutDir/RAxML_bipartitions.dendro 

# Kør Dendroscope og load .dendro fil

bash '/home/pato/dendroscope/Dendroscope'

# Åben RAxML_bipartitions.dendro 
# Dendroscope kan lave mange ændringer i udseende!



















# ANDET RAxML

# Laver mange træer
# This command will generate 20 ML trees on distinct starting trees and also print the tree with the best likelihood to a file called RAxML_bestTree.<Threads>. 
# -p er seed
raxmlHPC -m GTRGAMMA -p 12345 -s $DNA -# $NumberOfTrees -n $Threads


# Now we will want to get support values for this tree, so let's conduct a bootstrap search:
# Bootstrapper
# -b er seed for bootstrapping, For Rapid Bootstrapping, brug -x i stedet for -b
raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s dna.phy -n T12 
# Note that, RAxML also allows for automatically determining a sufficient number of bootstrap replicates, 
# in this case you would replace -# 100 by one of the bootstrap convergence criteria -# autoFC, -# autoMRE, -# autoMR, -# autoMRE_IGN.

# The nice thing about rapid bootstrapping is that it allows you to do a complete analysis (ML search + Bootstrapping) in one single step by typing
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s dna.phy -n T20
# If called like this RAxML will do 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values to you via one single program call.

# Bootstrap giver fil RAxML_bootstrap
# Now draw bipartitions on the best ML tree:
raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T12 -z RAxML_bootstrap.T12 -n T12

# This call will produce to output files that can be visualized with Dendroscope:  
# RAxML_bipartitions.T15 (support values assigned to nodes) and RAxML_bipartitionsBranchLabels.T15 
# (support values assigned to branches of the tree). Note that, for unrooted trees the correct representation 
# is actually the one with support values assigned to branches and not nodes of the tree!

# We can also use the Bootstrap replicates to build consensus trees, RAxML supports strict, majority rule, and extended majority rule consenus trees:
# 
#     strict consensus:            raxmlHPC -m GTRCAT -J STRICT -z RAxML_bootstrap.T14 -n T16
#     majority rule:                   raxmlHPC -m GTRCAT -J MR         -z RAxML_bootstrap.T14 -n T17
#     extended majority rule:  raxmlHPC -m GTRCAT -J MRE      -z RAxML_bootstrap.T14 -n T18


