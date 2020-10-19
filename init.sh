#!/bin/bash
set -e

mach=none

reinit=False
while getopts 'rm:' flag; do
    case "${flag}" in
    r) reinit=True ;;
    m) mach="$OPTARG" ;;
    *)
        echo "error"
        exit 1
        ;;
    esac
done

rootdir=$(pwd)

if [ "$reinit" == True ]; then
	if [ -f ._init_ ]; then
		rootdir1=$(sed -n 1p ._init_)
		sed -i "s|$rootdir1|_ROOTDIR_|g" scripts/*.sh
	else
		sed -i "s|$rootdir|_ROOTDIR_|g" scripts/*.sh
	fi
	
	echo _ROOTDIR_ > ._init_

	echo "Re-Initialized---"
else
	if [ -f ._init_ ]; then
		rootdir1=$(sed -n 1p ._init_)
		sed -i "s|$rootdir1|$rootdir|g" scripts/*.sh
	else
		sed -i "s|_ROOTDIR_|$rootdir|g" scripts/*.sh
	fi
	
	echo $rootdir > ._init_
	
	echo "RootDir Initialized---"
fi

if [ "$mach" != "none" ]; then
	echo $mach > $rootdir/bin/._machine_
	echo "Setting Machine as : $mach"
elif [ ! -z $rootdir/bin/._machine_ ]; then
	echo generic > $rootdir/bin/._machine_
	echo "Setting Machine as : generic"
fi

