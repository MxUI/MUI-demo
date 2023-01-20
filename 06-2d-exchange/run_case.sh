#!/bin/bash
rm -rf *.dat
mpirun -np 1 ./DomainCoarse : -np 1 ./DomainRefine
