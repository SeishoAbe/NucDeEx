#!/bin/bash

### FIXME ###
#TARGET_TMP="11B 11C 15N 15O"
#TARGET_TMP="15O 15N"
TARGET_TMP="11B"
ldmodel_tmp="1 2 3 4 5 6"
#############

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		mkdir fig/${TARGET}_ldmodel${ldmodel}
		mkdir output/${TARGET}_ldmodel${ldmodel}
		./bin/plot_decay $TARGET $ldmodel
	done
done
