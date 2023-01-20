#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3
HEADER_PATH	= -I/path/to/MUI

SCR = $(wildcard *.cpp)
EXE = $(SCR:.cpp=)

all: $(EXE)

% : %.cpp
	$(CC) $(CFLAGS) $(HEADER_PATH) $< -o $@

.PHONY: clean
	
clean:
	-rm -f $(EXE)