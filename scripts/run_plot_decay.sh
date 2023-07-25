#!/bin/bash

### FIXME ###
TARGET_TMP="11B 11C 15N 15O"
ldmodel_tmp="1 2 3"
parity_optmodall_tmp="0 1"
flag_jpi=1

#TARGET_TMP="14O 13O 12O 11O 14N 13N 12N 11N 10N 14C 13C 12C 11C 10C 9C 13B 12B 11B 10B 9B 8B 12Be 11Be 10Be 9Be 8Be 7Be 11Li 10Li 9Li 8Li 7Li 6Li"
#ldmodel_tmp=2
#parity_optmodall_tmp=1
#flag_jpi=0
#############

LogDir=log

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			FIGDIR=fig
			if [ $flag_jpi -eq 1 ] ; then
				if [ $TARGET = 11C ] || [ $TARGET = 11B ] ; then
					FIGDIR+=/12C
				elif [ $TARGET = 15N ] || [ $TARGET = 15O ] ; then
					FIGDIR+=/16O
				else 
					continue
				fi
			fi
			if [ ! -e $FIGDIR ] ; then
				mkdir $FIGDIR
			fi
			FIGDIR+=/${TARGET}_ldmodel${ldmodel}
			LogFile=$LogDir/log_plot_decay_${TARGET}_ldmodel${ldmodel}
			if [ $parity_optmodall -eq 1 ] ; then
				FIGDIR+="_parity_optmodall"
				LogFile+="_parity_optmodall"
			fi
			if [ $flag_jpi -eq 1 ] ; then 
				LogFile+="_flag_jpi"
			fi
			if [ ! -e $FIGDIR ] ; then
				mkdir $FIGDIR
			fi
			./bin/plot_decay $TARGET $ldmodel $parity_optmodall $flag_jpi > $LogFile 2>&1 &
		done
	done
done
