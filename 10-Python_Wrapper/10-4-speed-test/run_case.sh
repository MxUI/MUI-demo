#!/bin/bash

# Set the base directory to the directory where the script is located
cd "$(dirname "$0")"

Nt=10
Ni=100
Nj=100

START=$SECONDS

#NOTE: By using '-m mpi4py' mpi4py module runs the script handling properly
#      aborts. If a python process exit with failure all processes in MPI_COMM_WORLD       
#      will be aborted so hangs are avoided.

mpirun -np 1 python3 -m mpi4py fetcher.py $Nt $Ni $Nj :\
       -np 1 python3 -m mpi4py pusher.py $Nt $Ni $Nj 2>&1 | tee output.log

DURATION=$(( SECONDS - START ))

if (( $DURATION > 3600 )) ; then
    let "hours=DURATION/3600"
    let "minutes=(DURATION%3600)/60"
    let "seconds=(DURATION%3600)%60"
    echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)"
elif (( $DURATION > 60 )) ; then
    let "minutes=(DURATION%3600)/60"
    let "seconds=(DURATION%3600)%60"
    echo "Completed in $minutes minute(s) and $seconds second(s)"
elif (( $DURATION < 1 )) ; then
    echo "Completed in less than a second"
else
    echo "Completed in $DURATION seconds"
fi

#TODO: Segfaults with clang LLVM 9.1.0 randomly happened
#      Not sure if they are due to somebug in the binding code or the library itself.
#      I have not seen them with GCC 7.2.0.
