#!/bin/bash

queue="cccr-res"
WCLOCK="240"
JOBNAME="AMFI"
STDOUT="stdout"_$JOBNAME

submit_combine=True
nproc_combine=3
combine_only=True
combine_child_run=0
bypassrunning=True

EXE=_EXE_
RUNNCCP2R=_ROOTDIR_/exec/run_mppnccp2r/run_mppnccp2r


lsf=True
if ! [ -x "$(command -v bjobs)" ]; then
	lsf=False
fi

if [ ! "$bypassrunning" == "True" ];then
if [ "$lsf" == "True" ]; then
	alljobs=$(bjobs -noheader -o 'exec_cwd jobid job_name' 2>/dev/null | grep "$(pwd) ")
	if [ ! -z "$alljobs" ]; then
		echo "Some jobs are already running in this directory!"
	    echo "Please kill these jobs before submitting!"
		echo $alljobs
		exit 1
	fi
fi
fi

line1=$(sed -n '/&atmos_nml/,/\//p' input.nml | \
		sed -n '/layout/,/\//p' | \
		sed 's/=/ /g' | sed 's/,/ /g' | sed -e 's/\(.*\)/\L\1/')

found=0
npes=1
for strng in $line1; do
	if [ "$found" -gt "2" ]; then
		break
	elif [ "$found" -gt "0" ]; then
		found=$((found+1))
		npes=$((npes*strng))
	elif [[ "$strng" == layout ]]; then
		found=1
	fi
done

echo "NPES = "$npes 

if [ "$npes" -eq "1" ]; then
	submit_combine=False
fi

STDOUT=${STDOUT}

tfile=$(mktemp)
echo $tfile
cat <<EOF > $tfile

#!/bin/bash

#BSUB -q $queue
#BSUB -x
#BSUB -W ${WCLOCK}:00
#BSUB -J $JOBNAME
#BSUB -o $STDOUT
#BSUB -n $npes

export I_MPI_FABRICS=shm:dapl

ulimit -c unlimited
set -xu

mpirun -env OMP_NUM_THREADS 1 -n $npes $EXE
#mpirun -prepend-rank -env OMP_NUM_THREADS 1 -n $npes $EXE

EOF

rm -f $STDOUT

if [ "$submit_combine" == "True" ]; then
	if [ ! "$combine_only" == "True" ]; then
		echo "Removing any previous unprocessed output files."
		rm -f *.nc.????
	fi
fi

if [ "$combine_only" == "True" ]; then
	echo "Only combine"
	COND=""
else
  	cat $tfile
	if [ "$lsf" == "True" ]; then
		output=$(bsub < $tfile)
		echo $output
		jobid=$(echo $output | head -n1 | cut -d'<' -f2 | cut -d'>' -f1;)
		if [ "$jobid" -eq "$jobid" ] 2>/dev/null; then
		  	echo "Job submitted" 
		else
		  	exit 1
		fi
		COND="#BSUB -w started($jobid)"
	else
		bash $tfile
	fi
fi

if [ "$submit_combine" == "True" ]; then

WCLOCK=$((WCLOCK*2))

tfile=$(mktemp)
echo $tfile
cat <<EOF > $tfile

#!/bin/bash

#BSUB -q $queue
#BSUB -x
#BSUB -W ${WCLOCK}:00
#BSUB -J ${JOBNAME}_combine
#BSUB -o ${STDOUT}_combine
#$COND
#BSUB -n $nproc_combine

ulimit -c unlimited
set -xu

mpirun -n $nproc_combine $RUNNCCP2R \
<<< "&opts_nml removein=1, atmpes=$npes, child_run=$combine_child_run /"
EOF

rm -f ${STDOUT}_combine

cat $tfile
if [ "$lsf" == "True" ]; then
	bsub < $tfile
else 
	bash $tfile
fi

fi
