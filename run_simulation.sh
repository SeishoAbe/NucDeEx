#!/bin/bash

### FIXME ###
TARGET_TMP="11B 11C 15N 15O"
ldmodel_tmp="1 2 3"
parity_optmodall_tmp="0 1"
#TARGET_TMP="11B"
#ldmodel_tmp="1"
#parity_optmodall_tmp="0"
#############

LogDir=log

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			LogFile=$LogDir/log_${TARGET}_ldmodel${ldmodel}
			if [ $parity_optmodall -eq 1 ] ; then
				LogFile+="_parity_optmodall"
			fi
			./bin/simulation $TARGET $ldmodel $parity_optmodall > $LogFile 2>&1 &
		done
	done
done
