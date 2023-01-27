#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3

SCR = $(wildcard *.cpp)
EXE = $(SCR:.cpp=)

MUI_INC ?= 

default: $(EXE)

% : %.cpp
	$(CC) $(CFLAGS) -I$(MUI_INC) $< -o $@

clean:
	rm -f $(EXE) brownian.txt vortex.txt
