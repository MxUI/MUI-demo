#!/bin/bash

URI1="mpi://domain1/ifs"
URI2="mpi://domain2/ifs"

mpirun -np 1 hello ${URI1} 0.618 : -np 1 hello ${URI2} 1.414
