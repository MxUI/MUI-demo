 [MUI minimum required: (VERSION 2.0+)]
# MUI C Wrapper Demo: 3-pseudo-diffusion
This demo of the C wrapper of MUI Multiscale Universal Interface is a 3D scalar field pseudo diffusion case. The C++ `3D-pseudo-diffusion-fine` file generates two discontinuous domains (left domain and right domain viewed towards the negative z-axis direction) with a fine mesh density, while the C `3D-pseudo-diffusion-coarse` file generates one middle domain to "connect" the left and right domains with coarse mesh density. These three domains have the same dimension and align in the x-axis direction. The minimum x-axis boundary of the left domain has a scalar field source, which have a maximum value at the boundary centre and gradually decreases to outward. With the time increases, the scalar field diffuses along the positive x-axis direction. At each time instants, the left boundary of the middle domain fetches the field values from the right boundary of the left domain and the right boundary of the middle domain push field values to the left boundary of the right domain through MUI 2D interfaces with the RBF spatial sampler.

| Sketch of this demo |
|:------------------------------------------------------------------------------------------------------------------:|
| <img src='Resources/Sketch.jpg' width='1451px'/>                                                                           |

| Result at time-step = 200 |
|:------------------------------------------------------------------------------------------------------------------:|
|<img src='Resources/Result.jpg' width='1316px'/>|

This demo shows how to:
 1. Use the spatial sampler of Radial Basis Function between C++ and C codes;

This demo contains:
 1. 3D-pseudo-diffusion-fine.cpp & 3D-pseudo-diffusion-coarse.c: coupling domains; 
 2. CMakeLists.txt: CMake file;
 3. README.md

## Licensing

The source code is dual-licensed under both GPL v3 and Apache v2.

## Usage

To run the examples you need a `mpic++` wrapper with C++11 enabled backend, CMake V3.9+ installed, C wrapper of MUI compiled.
To run the demo:

```bash
bash run_case.sh
```
`paraview` will open automatically after running to visualise the results.

The integration of the scalar field:

The instantaneous integration of the scalar field of each boundary faces that perpendicular to the x-axis direction is in faceIntegrationD1_D3.txt and faceIntegrationD2.txt of the results folder.  

To clean the demo directory:

```bash
bash clean_case.sh
```

## Contact

Should you have any question, please do not hesitate to contact the MUI team at STFC
