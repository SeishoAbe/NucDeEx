#!/bin/bash

### FIXME ###
TARGET_TMP="11B 11C 15N 15O"
ldmodel_tmp="1 2 3 4 5 6"
parity_tmp="0 1"
optmodall_tmp="0 1"
#############

EXE=./bin/cout_input

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		for parity in $parity_tmp
		do
			for optmodall in $optmodall_tmp
			do
					$EXE $TARGET $ldmodel $parity $optmodall
			done
		done
	done
done
