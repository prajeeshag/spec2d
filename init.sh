#!/bin/bash
set -e

reinit=False
while getopts 'r' flag; do
    case "${flag}" in
    r) reinit=True ;;
    *)
        echo "error"
        exit 1
        ;;
    esac
done

rootdir=$(pwd)

if [ "$reinit" == True ]; then
	if [ -f ._init_ ]; then
		rootdir1=$(cat ._init_)
		sed -i "s|$rootdir1|_ROOTDIR_|g" scripts/*.sh
	else
		sed -i "s|$rootdir|_ROOTDIR_|g" scripts/*.sh
	fi
	
	echo _ROOTDIR_ > ._init_

	echo "Re-Initialized---"
else
	if [ -f ._init_ ]; then
		rootdir1=$(cat ._init_)
		sed -i "s|$rootdir1|$rootdir|g" scripts/*.sh
	else
		sed -i "s|_ROOTDIR_|$rootdir|g" scripts/*.sh
	fi
	
	echo $rootdir > ._init_
	
	echo "Initialized---"
fi
