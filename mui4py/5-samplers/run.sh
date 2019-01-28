#!/bin/bash
mpirun -np 1 python -m mpi4py speed_test.py dom1 1.0 : -np 1 python -m mpi4py speed_test.py dom2 2.0

