#!/bin/bash

TARGET=11C

cd /home/seisho/talys_1.96
. setup_talys_1.96.sh
cd -

OUTDIR=output_$TARGET

if [ ! -e $OUTDIR ] ; then
  echo "not output dir"
  exit
else 
  cd $OUTDIR
fi

FILE=energy.spin.parity
rm -v $FILE
ln -sv $TALYS_WORK_TABLES/energy_distribution/$FILE ./

FILE=energy
rm -v $FILE
ln -sv $TALYS_WORK_TABLES/energy_distribution/$FILE ./

talys < ../input_$TARGET > output_$TARGET
