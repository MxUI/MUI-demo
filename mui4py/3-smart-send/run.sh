#!/bin/bash
mpirun -np 1 python -m mpi4py pusher.py :\
       -np 4 python -m mpi4py fetcher.py 
