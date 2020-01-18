
machine=$(cat _ROOTDIR_/bin/._machine_)
source _ROOTDIR_/bin/env.$machine

alljobs=$(qstat -f  | grep "PBS_O_WORKDIR=$(pwd),")

workdir=$(echo $(pwd) | sed 's_/_\\/_g')

echo $alljobs
if [ ! -z "$alljobs" ]; then
	ids=$(qstat -f | grep -e "PBS_O_WORKDIR\|Job Id:" | sed -n "/${workdir},/{x;p;d;}; x" | awk '{print $3}')
	qstat -i $ids
	echo 
	echo
    echo "Some jobs are already running in this directory!"
    echo "Please kill these jobs before submitting!"
	echo 
	echo
	for id in $ids; do
		echo -n "Enter y and press [ENTER] to kill the Job id $id :"
		read inp
		if [ $inp == "y" ]; then
			qdel $id
		fi
	done

fi
