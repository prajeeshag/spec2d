#!/bin/bash

rootdir=$(pwd)

if [ ! -f ._init_ ]; then
	rootdir1=$(cat ._init_)
	sed -i "s|$rootdir1|$rootdir|g" scripts/*.sh
else
	sed -i "s|_ROOTDIR_|$rootdir|g" scripts/*.sh
fi

echo $rootdir > ._init_

echo "Initialized---"
