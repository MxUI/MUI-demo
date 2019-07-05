#!/bin/bash
#NOTE: By using '-m mpi4py' mpi4py module runs the script handling properly
#      aborts. If a python process exit with failure all processes in MPI_COMM_WORLD       
#      will be aborted so hangs are avoided.
mpirun -np 1 python3 -m mpi4py pusher.py dom1 :\
       -np 1 python3 -m mpi4py fetcher.py dom2
