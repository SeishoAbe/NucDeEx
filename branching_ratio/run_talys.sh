#!/bin/bash

### FIXME ###
TARGET_TMP="11B 11C 15N 15O"
FILE_TMP="energy.1.2.p.txt energy"
#############

OUTDIR=output
INDIR=input

cd /home/seisho/talys_1.96
. setup_talys_1.96.sh
cd -


if [ ! -e $OUTDIR ] ; then
  echo "not output dir"
  exit
else 
  cd $OUTDIR
fi

for FILE in $FILE_TMP
do
	ln -sv $TALYS_WORK_TABLES/energy_distribution/$FILE ./
done

for TARGET in $TARGET_TMP
do
	talys < ../$INDIR/input_$TARGET > output_$TARGET
done
