#!/bin/bash

prefix_tmp="SF_LFG_-14_1000_CCQE_C SF_LFG_-14_1000_CCQE_O SF_LFG_14_1000_CCQE_C SF_LFG_14_1000_CCQE_O  SF_LFG_-14_1000_NCQE_C SF_LFG_-14_1000_NCQE_O SF_LFG_14_1000_NCQE_C SF_LFG_14_1000_NCQE_O"

Exe=./bin/nuwro
LogDir=./log_nuwro

for prefix in $prefix_tmp
do
	LogFile=$LogDir/log_$prefix
	(time -f "CPU time %E: Memory %M KB" $Exe $prefix) > $LogFile 2>&1 &
done 
