#!/bin/bash
mpirun -np 2 python -m mpi4py code0.py dom0 :\
       -np 1 python -m mpi4py code1.py dom1 :\
       -np 2 python -m mpi4py code1.py dom2 :\
       -np 1 python -m mpi4py code2.py dom3 :\
       -np 1 python -m mpi4py code3.py dom4

