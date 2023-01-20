#!/bin/bash
mpirun -np 1 python3 -m mpi4py pusher.py :\
       -np 4 python3 -m mpi4py fetcher.py 
