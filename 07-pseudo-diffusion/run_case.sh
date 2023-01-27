#!/bin/bash

mpirun -np 1 3D-pseudo-diffusion-fine : \
       -np 1 3D-pseudo-diffusion-coarse
