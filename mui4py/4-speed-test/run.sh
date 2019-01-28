#!/bin/bash
Nt=10
Ni=100
Nj=100
mpirun -np 1 python -m mpi4py fetcher.py $Nt $Ni $Nj :\
       -np 1 python -m mpi4py pusher.py $Nt $Ni $Nj

#TODO: Segfaults with clang LLVM 9.1.0 randomly happened
#      Not sure if they are due to somebug in the binding code or the library itself.
#      I have not seen them with GCC 7.2.0.

