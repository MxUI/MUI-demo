#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -g -O3 -fopenmp

multidomain : multidomain.cpp
	$(CC) $(CFLAGS) $< -o $@
	@printf '\n\nTest example 1\n'
	@printf '    A\n   /|\\\n  i j k\n /  |  \\\nB   C   D\n'
	mpirun -np 1 ./$@ A i j k : -np 1 ./$@ B i : -np 1 ./$@ C j : -np 1 ./$@ D k
	@printf '\n\nTest example 2\n'
	@printf '    A*2\n   / \\\n  i   j\n /     \\\nB---k---C\n'
	mpirun -np 2 ./$@ A i j : -np 1 ./$@ B k i : -np 1 ./$@ C k j

clean:
	rm -f multidomain
