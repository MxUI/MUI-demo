#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3

SCR = $(wildcard *.cpp)
EXE = $(SCR:.cpp=)

MUI_INC ?= 
EIGEN_INC ?= 

default: $(EXE)

% : %.cpp
	$(CC) $(CFLAGS) -I$(MUI_INC) -I$(EIGEN_INC) -DUSE_RBF $< -o $@

clean:
	rm -rf $(EXE) coupling_results rbfCoarseMatrix rbfFineMatrix
