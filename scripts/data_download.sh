#!/bin/bash
set -e

datapath=https://www.dropbox.com/s/593vig1attxezbt/data_spec2d.tar.gz?dl=0 

if [ -d "/moes/home/prajeesh/spec2d" ]; then
	cd /moes/home/prajeesh/spec2d
	if [ -d "data" ]; then
		echo data directory already exist in /moes/home/prajeesh/spec2d
		exit
	fi

	if [ ! -f "data_spec2d.tar.gz" ]; then
		wget $datapath -O data_spec2d.tar.gz
	else
		echo "data tar file already present..!!"
	fi
	tar -zxvf data_spec2d.tar.gz
else
	echo /moes/home/prajeesh/spec2d does not exist, perhaps you might need to run init.sh first!
	exit 1
fi

