#!/bin/bash

FILE=energy
OUTDIR=output
INDIR=input
EXE_NUCLEUS=../read_nucleus.sh

cd /home/seisho/talys_1.96
. setup_talys_1.96.sh
cd -


if [ ! -e $OUTDIR ] ; then
  echo "not output dir"
  exit
else 
  cd $OUTDIR
fi

rm -v $FILE
ln -sv $TALYS_WORK_TABLES/energy_distribution/$FILE ./

$EXE_NUCLEUS | while read TARGET
do
	#echo $TARGET
	INPUTFILE=../$INDIR/input_$TARGET
	OUTPUTFILE=output_$TARGET

	if [ ! -e $INPUTFILE ] ; then
		echo "No inputfile: $INPUTFILE -> skip"
		continue
	fi
	echo "Process: $TARGET"
	talys < $INPUTFILE > $OUTPUTFILE
done
