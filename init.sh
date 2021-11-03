#!/bin/bash
set -e

mach='none'
avmach=''
for f in bin/env.*; do
	mc=$(echo $f | sed 's/bin\/env.//g')
	avmach="$avmach $mc"
done

usage(){
	echo ""
	echo "Usage: $0 -m machine_name"
	echo 
	echo "Available machines are: " $avmach
	exit
}

while getopts 'm:' flag; do
    case "${flag}" in
    m) mach="$OPTARG" ;;
    *) usage ;;
    esac
done

if [ -f .env ]; then
	echo "Error: .env file already exist!"
	exit 1
fi

rootdir=$(pwd)

if [[ "$mach" != "none" ]] && [[ "$avmach" == *"$mach"* ]]; then
	echo "export MACH=$mach" >> .env
	echo "Setting Machine as : $mach"
else
	usage
fi
	
echo "export rootdir=$rootdir" >> .env