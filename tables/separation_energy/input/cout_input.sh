#!/bin/bash


NUCLEUSTABLE=$TALYS_WORK_TABLES/nucleus/nucleus.txt

cat $NUCLEUSTABLE | while read line
do 
	check=${line:0:1}
	if [ $check = "#" ] ; then
		continue
	fi
	nucleus=`echo $line |cut -d ' ' -f 1`
	Z=`echo $line |cut -d ' ' -f 2`
	N=`echo $line |cut -d ' ' -f 3`
	#maxlevelsbin=`echo $line |cut -d ' ' -f 4`
	#echo $nucleus $Z $N $maxlevelsbin
	A=`expr $Z + $N`
	echo $nucleus $Z $N $A
	
	FILENAME=input_$nucleus
	{
		echo "projectile 0"
		echo "element $Z"
		echo "mass $A"
		echo "energy energy"
		echo "bins 0"
		echo "ejectiles g n p a d t h"
	} > $FILENAME


done
