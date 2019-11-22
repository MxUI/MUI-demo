# MUI Python Wrapper Demo: 5-Cpp-Python
This demo of the Python wrapper of MUI Multiscale Universal Interface is to exchange data between a C++ code and a Python code by using Python wrapper of MUI. 

This demo shows how to:
 1. Use the Python wrapper to exchange data between C++ and Python codes;

This demo contains:
 1. PUSHER_FETCHER_0.cpp & PUSHER_FETCHER_1.py: coupling codes; 
 2. Makefile;
 3. Allrun & Allclean: run control files;
 4. README.md

## Licensing

The source code is dual-licensed under both GPL v3 and Apache v2.

## Usage

To run the examples you need a `mpic++` wrapper with C++11 enabled backend, CMake V2.8+ installed, `Eigen3` installed, Python wrapper of MUI (mui4py) compiled.

To run the demo:

Complete the below line of the Makefile with the MUI code path

```bash
HEADER_PATH	= -I/path/to/MUI
```
and run the demo by

```bash
./Allrun
```

## Contact

Should you have any question, please do not hesitate to contact the MUI team at STFC
