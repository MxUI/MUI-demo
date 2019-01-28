#!/bin/bash
mpirun -np 1 python -m mpi4py data_types.py dom1  :\
       -np 1 python -m mpi4py data_types.py dom2
