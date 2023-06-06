#!/bin/bash

FILE=$TALYS_WORK_TABLES/nucleus/nucleus.txt

cat $FILE | while read line
do 
	check=${line:0:1}
	if [ $check = "#" ] ; then
		continue
	fi
	nucleus=${line%% *}
	echo $nucleus
done
