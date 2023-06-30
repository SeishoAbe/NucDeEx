#!/bin/bash

### FIXME ###
#TARGET_TMP="11B 11C 15N 15O"
#ldmodel_tmp="1 2 3"
#parity_optmodall_tmp="0 1"
TARGET_TMP="15N"
ldmodel_tmp="1"
parity_optmodall_tmp="0"
#############
FILE_TMP="energy.1.2.p.txt energy"

#ldmodel_tmp="1 2 3 4 5 6"
# 4-6 (microscopic is not available)

INDIR=input

cd /home/seisho/talys_1.96
. setup_talys_1.96.sh
cd -

for TARGET in $TARGET_TMP
do
	OUTDIR=output/$TARGET
	if [ ! -e $OUTDIR ] ; then
		mkdir $OUTDIR
	fi
	cd $OUTDIR

	for FILE in $FILE_TMP
	do
		ln -sv $TALYS_WORK_TABLES/energy_distribution/$FILE ./
	done

	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			INFILE=$TALYS_WORK/$INDIR/input_${TARGET}_ldmodel${ldmodel}
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
