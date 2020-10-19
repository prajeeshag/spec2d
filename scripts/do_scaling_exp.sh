#!/bin/bash

set -e

NLAT=94
TRUNC=62
minpes=8
maxpes=1000

#mv exp/valid_pe_layouts_${NLAT}_T${TRUNC} valid_pe_layouts

#mv exp/input.nml .

#mv exp/run_amfi_*.sh run_amfi.sh

layout_list=$(awk -v pes1=$minpes -v pes2=$maxpes '{ if ($2 == "pes:-" && $1 >= pes1 && $1 <= pes2) print $0}' valid_pe_layouts | \
	awk '{ $1="";$2=""; print $0}' | sed -e 's/;/\n/g' -e 's/\[//g' -e 's/\]//g' | awk 'NF')


for layout in  $layout_list; do
	xpe=$(echo $layout | awk -F "," '{print $1}')
	ype=$(echo $layout | awk -F "," '{print $2}')
	npes=$((xpe*ype))
	jobnm=AMFI_${NLAT}_${TRUNC}_${npes}_${xpe}x${ype}
	echo $jobnm
	rm -rf $jobnm
    mkdir $jobnm
	ln -s $(pwd)/exp/* $jobnm/
	rm $jobnm/RESTART
	mkdir $jobnm/RESTART
	sed "s/_LAYOUT_/$layout/" input.nml > $jobnm/input.nml
	sed "s/_JOBNAME_/$jobnm/" run_amfi.sh > $jobnm/run_amfi.sh
	cd $jobnm
	bash ./run_amfi.sh
	cd -
done
