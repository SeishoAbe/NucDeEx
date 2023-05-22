#!/bin/bash

cd /home/seisho/talys_1.96
. setup_talys_1.96.sh
cd -

if [ ! -e output ] ; then
  echo "not output dir"
  exit
else 
  cd output
fi

FILE=energy.spin.parity
rm -v $FILE
ln -sv ../../../tables/energy_distribution/$FILE ./

FILE=energy
rm -v $FILE
ln -sv ../../../tables/energy_distribution/$FILE ./

talys < ../input > output
