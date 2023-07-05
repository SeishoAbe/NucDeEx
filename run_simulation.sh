#!/bin/bash

### FIXME ###
TARGET_TMP="11B 11C 15N 15O"
ldmodel_tmp="1 2 3"
parity_optmodall_tmp="0 1"
flag_jpi=1

#TARGET_TMP="11B 11C 15N 15O"
#ldmodel_tmp="2"
#parity_optmodall_tmp="1"
#############

LogDir=log

for TARGET in $TARGET_TMP
do
	FIGDIR=fig_sim
	OUTDIR=sim_out
	if [ $flag_jpi -eq 1 ] ; then
		if [ $TARGET = 11C ] || [ $TARGET = 11B ] ; then
			FIGDIR+=/12C
			OUTDIR+=/12C
		elif [ $TARGET = 15N ] || [ $TARGET = 15O ] ; then
			FIGDIR+=/16O
			OUTDIR+=/16O
		else 
			continue
		fi
	fi
	if [ ! -e $FIGDIR ] ; then
		mkdir $FIGDIR
	fi
	if [ ! -e $OUTDIR ] ; then
		mkdir $OUTDIR
	fi
	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			LogFile=$LogDir/log_simulation_${TARGET}_ldmodel${ldmodel}
			if [ $parity_optmodall -eq 1 ] ; then
				LogFile+="_parity_optmodall"
			fi
			(time -f "Memory %M KB" ./bin/simulation $TARGET $ldmodel $parity_optmodall $flag_jpi) > $LogFile 2>&1 &
		done
	done
done
