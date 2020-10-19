#!/bin/bash

set -e

NLAT=94
TRUNC=62
minpes=8
maxpes=1000

layout_list=$(awk -v pes1=$minpes -v pes2=$maxpes '{ if ($2 == "pes:-" && $1 >= pes1 && $1 <= pes2) print $0}' valid_pe_layouts | \
	awk '{ $1="";$2=""; print $0}' | sed -e 's/;/\n/g' -e 's/\[//g' -e 's/\]//g' | awk 'NF')

for layout in  $layout_list; do
	xpe=$(echo $layout | awk -F "," '{print $1}')
	ype=$(echo $layout | awk -F "," '{print $2}')
	npes=$((xpe*ype))
	jobnm=AMFI_${NLAT}_${TRUNC}_${npes}_${xpe}x${ype}
	echo $jobnm
	cd $jobnm
	bash ./run_amfi.sh
	cd -
done
