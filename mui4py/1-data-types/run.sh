#!/bin/bash
mpirun -np 1 python3 -m mpi4py data_types.py dom1  :\
       -np 1 python3 -m mpi4py data_types.py dom2
