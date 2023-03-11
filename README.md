# MUI Demo
This repository contains demos and example solvers that make use of the MUI Multiscale Universal Interface to perform multiscale/multiphysics computations. The MUI primary repository can be found at [https://github.com/MxUI/MUI](https://github.com/MxUI/MUI).

## Licensing

The source code is dual-licensed under both GPL v3 and Apache v2.

## Run Demos

To run the demos you need a copy of the [MUI code](https://github.com/MxUI/MUI), a `mpic++` wrapper with C++11 enabled backend. You also need to install C, FORTRAN and/or Python wrappers of MUI if you want to run demos for wrappers.

To run demos:

```bash
cd <repository_directory>/<demo_directory>
./run_case.sh
```
A log file "output.log" will be generated.

In case you encounter an error related to `mui.h` not being found during compilation with CMake, you can resolve this issue by specifying the path to MUI using the command:

```bash
./run_case.sh /path/to/MUI
```

For demo 02-heat, you can view the time evolution of the simulation by running the vizmulti.m MATLAB script and dragging the slider bar in the figure.

```bash
matlab -r Resources/vizmulti
```
For demos 03-heat-sph-fdm and 04-bd-ns, you can view the result by using gnuplot with the following command. Also checkout flow.png in demo 04-bd-ns for a visualisation of the flow field:

```bash
gnuplot -p Resources/plot.gp
```

For demos 09-algorithms, Paraview state file will be automatically executed to display the results after running the code.

To clean the demo directory:

```bash
cd <repository_directory>/<demo_directory>
./clean_case.sh
```

## Publication

**Tang** Y.-H., **Kudo**, S., **Bian**, X., **Li**, Z., & **Karniadakis**, G. E. Multiscale Universal Interface: A Concurrent Framework for Coupling Heterogeneous Solvers, *Journal of Computational Physics*, **2015**, 297.15, 13-31.

## Contact

Should you have any question please do not hesitate to contact the developers, a list can be found within the MxUI GitHub organisation pages [https://github.com/MxUI/MUI](https://github.com/MxUI/MUI)
