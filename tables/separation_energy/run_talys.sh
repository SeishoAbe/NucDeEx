#!/bin/bash

FILE=energy
OUTDIR=output
INDIR=input
NUCLEUSTABLE=$TALYS_WORK_TABLES/nucleus/nucleus.txt

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

cat $NUCLEUSTABLE | while read line
do
	check=${line:0:1}
	if [ $check = "#" ] ; then
		continue
	fi
	TARGET=`echo $line |cut -d ' ' -f 1`
	INPUTFILE=../$INDIR/input_$TARGET
	OUTPUTFILE=output_$TARGET

	if [ ! -e $INPUTFILE ] ; then
		echo "No inputfile: $INPUTFILE -> skip"
		continue
	fi
	echo "Process: $TARGET"
	talys < $INPUTFILE > $OUTPUTFILE
done
