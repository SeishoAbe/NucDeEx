#!/bin/bash

### FIXME ###
TARGET_TMP="11B 11C 15N 15O"
#############

EXE=./bin/cout_input

for TARGET in $TARGET_TMP
do
	$EXE $TARGET
done
