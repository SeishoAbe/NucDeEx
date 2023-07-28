#! /bin/bash

### FIXME ###
target_tmp="C O"
tune_tmp="G18_10a_02_11a G18_10a_02_11b G18_10b_02_11a G18_10b_02_11b G23_10b_02_11a G23_10b_02_11b"
#tune_tmp="G18_10a_02_11a G18_10a_02_11b"
#tune_tmp="G18_10b_02_11a G18_10b_02_11b"
#tune_tmp="G23_10b_02_11a G23_10b_02_11b"

List_tmp="NCEL CCQE"
#############

Exe=./bin/genie
LogDir=./log_genie

for List in $List_tmp
do
	for pdg in 14 -14 
	do
		for target in $target_tmp
		do
			for tune in $tune_tmp
			do
				LogFile=$LogDir/log_${List}_${pdg}_${target}_${tune}
			(time -f "CPU time %E: Memory %M KB"	$Exe $List $pdg $target $tune) > $LogFile 2>&1 &
			done 
		done 
	done 
done 
