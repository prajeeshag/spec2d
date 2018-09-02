#!/bin/bash

queue="cccr-res"
WCLOCK="240:00"
JOBNAME="AMFI"
STDOUT="stdout"
EXE=/moes/home/prajeesh/spec2d/exec/spec2d/spec2d.exe


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

tfile=$(mktemp)

STDOUT=${STDOUT}_${npes}

rm -f $STDOUT

cat <<EOF > $tfile

#!/bin/bash

#BSUB -q $queue
#BSUB -x
#BSUB -W $WCLOCK
#BSUB -J $JOBNAME
#BSUB -o $STDOUT
#BSUB -n $npes

export I_MPI_FABRICS=shm:dapl

ulimit -c unlimited
set -xu

mpirun -env OMP_NUM_THREADS 1 -n $npes $EXE
#mpirun -prepend-rank -env OMP_NUM_THREADS 1 -n $npes $EXE

EOF

bsub < $tfile
