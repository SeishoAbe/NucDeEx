#!/bin/bash

FILE=energy
OUTDIR=output
INDIR=input
NUCLEUSTABLE=$NUCDEEX_TABLES/nucleus/nucleus.txt

if [ ! -e $OUTDIR ] ; then
  echo "not output dir"
  exit
else 
  cd $OUTDIR
fi

rm -v $FILE
ln -sv $NUCDEEX_TABLES/energy_distribution/$FILE ./

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
