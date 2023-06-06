#!/bin/bash

EXE=./cout_separation_energy
EXE_NUCLEUS=./read_nucleus.sh
OUTDIR=output

$EXE_NUCLEUS | while read TARGET
do
	#check file
	OUTPUTFILE=output/output_$TARGET

	if [ ! -e $OUTPUTFILE ] ; then
		echo "No ouptut file: $OUTPUTFILE -> skip"
		continue
	fi
	echo "Process: $TARGET"
	$EXE $TARGET
done
