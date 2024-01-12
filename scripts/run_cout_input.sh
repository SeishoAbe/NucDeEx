#!/bin/bash

### FIXME ###

### single nucleon hole
#TARGET_TMP="11B 11C 15N 15O"
#ldmodel_tmp="1 2 3" #4 5 6
#parity_optmodall_tmp="0 1"
#flag_jpi=1

### multi nucleon hole
#TARGET_TMP="14O 13O 12O 11O 14N 13N 12N 11N 10N 14C 13C 12C 11C 10C 9C 13B 12B 11B 10B 9B 8B 12Be 11Be 10Be 9Be 8Be 7Be 11Li 10Li 9Li 8Li 7Li 6Li"
TARGET_TMP="16O 16N 16F"
ldmodel_tmp=2
parity_optmodall_tmp=1
flag_jpi=0
#############

EXE=./bin/cout_input

for TARGET in $TARGET_TMP
do
	for ldmodel in $ldmodel_tmp
	do
		for parity_optmodall in $parity_optmodall_tmp
		do
			$EXE $TARGET $ldmodel $parity_optmodall $parity_optmodall $flag_jpi
		done
	done
done
