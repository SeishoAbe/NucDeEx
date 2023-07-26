#!/bin/bash

prefix_tmp="neut_1GeV_numu_CCQE_C_MDLQE422 neut_1GeV_numub_CCQE_C_MDLQE422 neut_1GeV_numu_CCQE_O_MDLQE422_NUCDEXITE0 neut_1GeV_numub_CCQE_O_MDLQE422_NUCDEXITE0"

Exe=./bin/neut
LogDir=./log_neut

for prefix in $prefix_tmp
do
	LogFile=$LogDir/log_$prefix
	(time -f "CPU time %E: Memory %M KB" $Exe $prefix) > $LogFile 2>&1 &
done 
