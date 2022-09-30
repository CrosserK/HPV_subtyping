#!/bin/bash

# Rydder midlertidige filer op efter endt FindHPVMutations.sh kørsel

# Henter options
Options=$(grep -v "^-" ../../FindHPVMutations_options.txt)

Name=$(grep "RunName=" <<< $Options | sed 's/RunName=//g')
MainF=$(grep "MainF=" <<< $Options | sed 's/MainF=//g')

if [ ${#Name} -gt 0 ] && [ ${#MainF} -gt 0 ]; then
rm -I $MainF/FASTQfiles_${Name}*.txt
rm -I --recursive -f $MainF/GenotypeCalls/${Name}*
rm -I --recursive $MainF/VirStrain_run/${Name}*
fi


