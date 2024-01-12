#!/bin/bash

### FIXME ###
#TARGET_TMP="11B 11C 15N 15O"
#ldmodel_tmp="1 2 3"
#parity_optmodall_tmp="0 1"
#flag_jpi=1

#TARGET_TMP="14O 13O 12O 11O 14N 13N 12N 11N 10N 14C 13C 12C 11C 10C 9C 13B 12B 11B 10B 9B 8B 12Be 11Be 10Be 9Be 8Be 7Be 11Li 10Li 9Li 8Li 7Li 6Li"
TARGET_TMP="16O"
ldmodel_tmp=2
parity_optmodall_tmp=1
flag_jpi=0
#############

FILE_TMP="energy.1.2.p energy"
#ldmodel_tmp="1 2 3 4 5 6" # 4-6 (microscopic is not available)

INDIR_PREFIX=input
OUTDIR_PREFIX=output

for TARGET in $TARGET_TMP
do
	INDIR=$INDIR_PREFIX
	OUTDIR=$OUTDIR_PREFIX
	if [ $flag_jpi -eq 1 ] ; then
		if [ $TARGET = 11C ] || [ $TARGET = 11B ] ; then
			INDIR+=/12C
			OUTDIR+=/12C
		elif [ $TARGET = 15N ] || [ $TARGET = 15O ] ; then
			INDIR+=/16O
			OUTDIR+=/16O
		else 
			continue
		fi
	fi

	if [ ! -e $OUTDIR ] ; then
		mkdir $OUTDIR
	fi
	cd $OUTDIR

	for FILE in $FILE_TMP
	do
		ln -sv $NUCDEEX_TABLES/energy_distribution/$FILE ./
	done

	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			INFILE=$NUCDEEX_ROOT/$INDIR/input_${TARGET}_ldmodel${ldmodel}
			OUTFILE=output_${TARGET}_ldmodel${ldmodel}
			if [ $parity_optmodall -eq 1 ] ; then
				INFILE+="_parity_optmodall"
				OUTFILE+="_parity_optmodall"
			fi
			talys < $INFILE > $OUTFILE
		done
	done
	cd -
done
