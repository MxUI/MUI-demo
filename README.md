# MUI Demo
This repository contains demos and example solvers that make use of the MUI Multiscale Universal Interface to perform multiscale/multiphysics computations. The MUI primary repository can be found at [https://github.com/yhtang/MUI](https://github.com/yhtang/MUI).

## Licensing

The source code is dual-licensed under both GPL v3 and Apache v2.

## Installation

To run the examples you need a `mpic++` wrapper with C++11 enabled backend. In
order to run the first example run

```bash
# clone the repository
cd <repository_directory>
git submodule update --init  # to obtain MUI
cd 0-hello-world
make
mpirun -np 1 ./hello mpi://domain1/ifs 0.618 : -np 1 ./hello mpi://domain2/ifs 1.414
```

## Publication

**Tang** Y.-H., **Kudo**, S., **Bian**, X., **Li**, Z., & **Karniadakis**, G. E. Multiscale Universal Interface: A Concurrent Framework for Coupling Heterogeneous Solvers, *Journal of Computational Physics*, **2015**, 297.15, 13-31.

## Contact

Should you have any question please do not hesitate to contact yuhang_tang at brown dot edu
