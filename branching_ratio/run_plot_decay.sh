#!/bin/bash

### FIXME ###
#TARGET=11B
TARGET=11C
#############

mkdir fig/$TARGET

./bin/plot_decay $TARGET
