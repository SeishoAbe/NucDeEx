#!/bin/bash

EXE=./bin/cout_separation_energy
OUTDIR=output
NUCLEUSTABLE=$NUCDEEX_TABLES/nucleus/nucleus.txt

cat $NUCLEUSTABLE | while read line
do
	check=${line:0:1}
	if [ $check = "#" ] ; then
		continue
	fi
	TARGET=`echo $line |cut -d ' ' -f 1`
	#check file
	OUTPUTFILE=output/output_$TARGET

	if [ ! -e $OUTPUTFILE ] ; then
		echo "No ouptut file: $OUTPUTFILE -> skip"
		continue
	fi
	echo "Process: $TARGET"
	$EXE $TARGET
done
