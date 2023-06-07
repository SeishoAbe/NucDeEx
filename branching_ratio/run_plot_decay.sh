#!/bin/bash

### FIXME ###
#TARGET_TMP="11B 11C 15N 15O"
#TARGET_TMP="15O"
TARGET_TMP="11B"
#############

for TARGET in $TARGET_TMP
do
	mkdir fig/$TARGET
	mkdir output/$TARGET
	./bin/plot_decay $TARGET
done
