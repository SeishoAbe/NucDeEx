#!/bin/bash

### FIXME ###
TARGET_TMP="11B"
ldmodel_tmp="1 2 3 4 5 6"
#############

LogDir=log

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		LogFile=$LogDir/log_${TARGET}_ldmodel${ldmodel}
		./bin/simulation $TARGET $ldmodel > $LogFile 2>&1 &
	done
done
