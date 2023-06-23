#!/bin/bash

### FIXME ###
#TARGET_TMP="11B 11C 15N 15O"
TARGET_TMP="15O 15N"
#TARGET_TMP="11B 11C"
ldmodel_tmp="1 2 3"
parity_optmodall_tmp="0 1"
#############

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			FIGDIR=fig/${TARGET}_ldmodel${ldmodel}
			if [ $parity_optmodall -eq 1 ] ; then
				FIGDIR+="_parity_optmodall"
			fi
			mkdir $FIGDIR
			./bin/plot_decay $TARGET $ldmodel $parity_optmodall
		done
	done
done
