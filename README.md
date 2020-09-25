# MUI Demo
This repository contains demos and example solvers that make use of the MUI Multiscale Universal Interface to perform multiscale/multiphysics computations. The MUI primary repository can be found at [https://github.com/MxUI/MUI](https://github.com/MxUI/MUI).

## Licensing

The source code is dual-licensed under both GPL v3 and Apache v2.

## Installation

To run the examples you need a `mpic++` wrapper with C++11 enabled backend.

To run the first example 0-hello-world:

```bash
# clone the repository
cd <repository_directory>
git submodule update --init  # to obtain MUI
cd 0-hello-world
make
mpirun -np 1 ./hello mpi://domain1/ifs 0.618 : -np 1 ./hello mpi://domain2/ifs 1.414
```

To run 2-heat:

```bash
cd <repository_directory>/1-heat
make
mpirun -np 1 ./heat-coarse : -np 1 ./heat-fine
matlab -r vizmulti
```

You can view the time evolution of the simulation by running the vizmulti.m MATLAB script and dragging the slider bar in the figure.

To run 4-bd-ns:

```bash
cd <repository_directory>/4-bd-ns
make
mpirun -np 1 ./brownian : -np 4 ./vortex
gnuplot -p plot.gp
```

Also checkout flow.png for a visualization of the flow field.

To run 5-multi-domain:

```bash
cd <repository_directory>/5-multi-domain
make
```

## Publication

**Tang** Y.-H., **Kudo**, S., **Bian**, X., **Li**, Z., & **Karniadakis**, G. E. Multiscale Universal Interface: A Concurrent Framework for Coupling Heterogeneous Solvers, *Journal of Computational Physics*, **2015**, 297.15, 13-31.

## Contact

Should you have any question please do not hesitate to contact the developers, a list can be found within the MxUI GitHub organisation pages [https://github.com/MxUI/MUI](https://github.com/MxUI/MUI)
