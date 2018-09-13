#!/bin/bash

queue="cccr-res"
WCLOCK="240"
JOBNAME="AMFI"
STDOUT="stdout"_$JOBNAME

# submit_combine - if "False" do not submit postprocessing
submit_combine=True
# nproc_combine - number of processors for postproc
nproc_combine=8
# combine_only - if "True" submit postprocessing only
combine_only=False
# combine_child_run - if "1" postprocessing is submitted as child run of model run
combine_child_run=1

EXE=_EXE_
RUNNCCP2R=_ROOTDIR_/exec/run_mppnccp2r/run_mppnccp2r












#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------- 


if [ ! "$submit_combine" == "True" ]; then
    combine_only=False
fi

lsf=True
if ! [ -x "$(command -v bjobs)" ]; then
	lsf=False
fi

bypassrunning=False
if [ "$combine_only" == "True" ]; then
	if [ "$combine_child_run" -eq "1" ]; then
		bypassrunning=True
	fi
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
	rm -f $STDOUT
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
		COND="BSUB -w started($jobid)"
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
