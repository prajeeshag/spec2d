#!/bin/bash

queue="cccr-res"
WCLOCK="240"
JOBNAME="AMFI"
STDOUT="stdout"_$JOBNAME
submit_combine=True


EXE=/moes/home/prajeesh/spec2d/exec/spec2d/spec2d.exe
RUNNCCP2R=/moes/home/prajeesh/spec2d/exec/run_mppnccp2r/run_mppnccp2r


alljobs=$(bjobs -noheader -o 'exec_cwd jobid job_name' 2>/dev/null | grep "$(pwd) ")
if [ ! -z "$alljobs" ]; then
	echo "Some jobs are already running in this directory!"
    echo "Please kill these jobs before submitting!"
	echo $alljobs
	exit 1
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

STDOUT=${STDOUT}_${npes}

tfile=$(mktemp)
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
	echo "Removing any previous unprocessed output files."
	rm -f *.nc.????
fi

output=$(bsub < $tfile)

echo $output

jobid=$(echo $output | head -n1 | cut -d'<' -f2 | cut -d'>' -f1;)

if [ "$jobid" -eq "$jobid" ] 2>/dev/null; then
  	echo "Job submitted" 
else
  	exit 1
fi


if [ "$submit_combine" == "True" ]; then

WCLOCK=$((WCLOCK*2))

tfile=$(mktemp)
cat <<EOF > $tfile

#!/bin/bash

#BSUB -q $queue
#BSUB -x
#BSUB -W ${WCLOCK}:00
#BSUB -J ${JOBNAME}_combine
#BSUB -o ${STDOUT}_combine
#BSUB -w started($jobid)
#BSUB -n 1

ulimit -c unlimited
set -xu

mpirun -n 1 $RUNNCCP2R <<< "&opts_nml removein=T, atmpes=$npes, child_run=T /"

EOF

rm -f ${STDOUT}_combine

bsub < $tfile

fi
