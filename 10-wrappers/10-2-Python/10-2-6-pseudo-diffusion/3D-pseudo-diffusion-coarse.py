"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library                       #
#                                                                            #
# Copyright (C) 2023 W. Liu                                                  #
#                                                                            #
# This software is jointly licensed under the Apache License, Version 2.0    #
# and the GNU General Public License version 3, you may use it according     #
# to either.                                                                 #
#                                                                            #
# ** Apache License, version 2.0 **                                          #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# You may obtain a copy of the License at                                    #
#                                                                            #
# http://www.apache.org/licenses/LICENSE-2.0                                 #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
# ** GNU General Public License, version 3 **                                #
#                                                                            #
# This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by       #
# the Free Software Foundation, either version 3 of the License, or          #
# (at your option) any later version.                                        #
#                                                                            #
# This program is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
# GNU General Public License for more details.                               #
#                                                                            #
# You should have received a copy of the GNU General Public License          #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
##############################################################################
#
# @file 3D-pseudo-diffusion-coarse.py
# @author W. Liu
# @date 04 April 2023
# @brief Coarse (middle) domain of the 3D pseudo diffusion case
#
"""

from mpi4py import MPI
import mpi4py
import datetime
import numpy as np
import time
import os

# Include MUI header file and configure file 
import mui4py

# Common world claims 
MUI_COMM_WORLD = mui4py.mpi_split_by_app()

# Declare MPI ranks
rank = MUI_COMM_WORLD.Get_rank()

# Declare MPI size
size = MUI_COMM_WORLD.Get_size()

# Define MUI dimension
dimensionMUI = 2

# Define the name of push/fetch values
name_pushA = "coarseFieldA"
name_fetchA  = "fineFieldA"

# Define MUI push/fetch data types
data_types_push = {name_pushA: mui4py.FLOAT64}
data_types_fetch = {name_fetchA: mui4py.FLOAT64} 

# MUI interface creation
domain = "coarseDomain"
config2d = mui4py.Config(dimensionMUI, mui4py.FLOAT64)
iface = ["interface2D01", "interface2D02"]
MUI_Interfaces = mui4py.create_unifaces(domain, iface, config2d) 
MUI_Interfaces["interface2D01"].set_data_types(data_types_fetch)
MUI_Interfaces["interface2D02"].set_data_types(data_types_push)

print("mpi4py.get_config(): ", mpi4py.get_config(), "\n")

# Define the forget steps of MUI to reduce the memory
forgetSteps = int(5)
# Define the synchronised boolen for MUI smart send
synchronised=False

#Define parameters of the RBF sampler
rSampler = 0.8                                  # Define the search radius of the RBF sampler (radius size should be balanced to try and maintain)
basisFunc = 1                                   # Specify RBF basis function 0-Gaussian; 1-WendlandC0; 2-WendlandC2; 3-WendlandC4; 4-WendlandC6
conservative = True                             # Enable conservative OR consistent RBF form
cutOff = 1e-9                                   # Cut-off value for Gaussian RBF basis function
smoothFunc = False                              # Enable/disable RBF smoothing function during matrix creation
writeMatrix = True                              # Enable/disable writing of the matrix (if not reading)
cgSolveTol = 1e-6;                              # Conjugate Gradient solver tolerance
cgMaxIter = 500;                                # Conjugate Gradient solver maximum iterations (-1 = value determined by tolerance)
preconditioner = 1;                             # Preconditioner of Conjugate Gradient solver
pouSize = 50;                                   # RBF Partition of Unity patch size
rbfMatrixFolder = "rbfCoarseMatrix" + str(rank) # Output folder for the RBF matrix files

# Create the RBF Matrix Folder includes all intermediate folders if don't exists
if not os.path.exists(rbfMatrixFolder):
    os.makedirs(rbfMatrixFolder)
    print("Folder " , rbfMatrixFolder, " created ")
else:
    print("Folder " , rbfMatrixFolder,  " already exists")

# Setup diffusion rate
dr = 0.5

# Setup time steps
steps = int(200)

# Setup output interval
outputInterval  = int(1)

# Domain discretization
# total number of grid points in x axis
Ntx = int(9)
# total number of grid points in y axis
Nty = int(9)
# total number of grid points in z axis
Ntz = int(9)

# number of grid points in x axis per MPI rank
Nx = Ntx
# number of grid points in y axis per MPI rank
Ny = Nty
# number of grid points in z axis per MPI rank
Nz = int(0)

if (rank < (Ntz % size)):
    Nz = int(Ntz/size + 1)
else:
    Nz = int(Ntz/size)

# total number of points
Npoints = Nx*Ny*Nz
assert Npoints > 0

# Geometry info
# length (x-axis direction) of the geometry
lx = 1.0
# length (y-axis direction) of the geometry
ly = 1.0
# length (z-axis direction) of the geometry
lzpr = 0.0

if(size == 2):
    if((Ntz%2)==0):
        lpz = 1.0/(float(Ntz)-1.0)
        lpz_half = lpz/2.0
        lzpr = (1.0 / float(size))-lpz_half
    else:
        if(rank == 0):
            lzpr = 1.0 / float(size)
        else:
            lpz = 1.0/(float(Ntz)-1.0)
            lzpr = (1.0 / float(size))-lpz
elif (size == 1):
    lzpr = 1.0

# origin coordinate (x-axis direction) of the geometry
x0 = 1.0
# origin coordinate (y-axis direction) of the geometry
y0 = 0.0
# origin coordinate (z-axis direction) of the geometry
z0pr = 0.0

if(rank == 0):
    z0pr = 0.0
elif (rank == 1):
    if((Ntz%2)==0):
        lpz = 1.0/(float(Ntz)-1.0)
        lpz_half = lpz/2.0
        z0pr = (1.0 / float(size))+lpz_half
    else:
        lpz = 1.0/(float(Ntz)-1.0)
        z0pr = (1.0 / float(size))+lpz

# Tolerance of data points
tolerance = (lx/(Nx-1))*0.5

# Store point coordinates
points = np.zeros((Npoints, 3))
c_0 = 0
for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            points[c_0] = [(x0+(lx/(Nx-1))*i), (y0+(ly/(Ny-1))*j), (z0pr+(lzpr/(Nz-1))*k)]
            c_0 += 1

# Generate initial pseudo scalar field
scalar_field_A = np.zeros(Npoints)
for i in range(Npoints):
    scalar_field_A[i] = 0.0

# Declare list to store mui::point2d
point2dList = []

# Store mui::point2d that located in the fetch interface
c_0 = 0
for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            if ((points[c_0][0]-x0) <= tolerance):
                point_fetch = MUI_Interfaces["interface2D01"].Point([points[c_0][1], points[c_0][2]])
                point2dList.append(point_fetch)
            c_0 += 1

# Define and announce MUI send/receive span
send_span = mui4py.geometry.Box([y0, z0pr], [(y0+ly), (z0pr+lzpr)])
MUI_Interfaces["interface2D02"].announce_send_span(0, (steps+1), send_span, synchronised)

# Spatial/temporal samplers
t_sampler = mui4py.TemporalSamplerExact()
s_sampler = mui4py.SamplerRbf(rSampler, point2dList, basisFunc, conservative, smoothFunc, writeMatrix, rbfMatrixFolder, cutOff, cgSolveTol, cgMaxIter, pouSize, preconditioner, MUI_COMM_WORLD)

# Commit ZERO step
MUI_Interfaces["interface2D02"].commit(0)

# Output the initial pseudo scalar field
outputFileName = "coupling_results" + str(rank) + "/scalar_field_coarse_0.csv"
outputFile = open(outputFileName,"w+")
outputFile.write("\"X\",\"Y\",\"Z\",\"scalar_field_A\"\n")
c_0 = 0
for k in range(Nz):
    for j in range(Ny):
        for i in range(Nx):
            outputFile.write("%f,%f,%f,%f\n" % (points[c_0][0],points[c_0][1],points[c_0][2],scalar_field_A[c_0]))
            c_0 += 1
outputFile.close() 

# Declare integrations of the pseudo scalar field at boundaries
# Integration of the pseudo scalar field at the left boundary of Domain 2
intFaceLD2 = 0.0
# Integration of the pseudo scalar field at the right boundary of Domain 2
intFaceRD2 = 0.0

# Define output files for boundary integrations
outputIntegrationName = "coupling_results" + str(rank) + "/faceIntegrationD2.txt"
outputIntegration = open(outputIntegrationName,"w+")
outputIntegration.write("\"t\",\"intFaceLD2\",\"intFaceRD2\"\n")
outputIntegration.close()

# Begin time loops
for t in range(1, (steps+1)):

    print("\n")
    print("{Coarse Domain} ", t," Step \n")

    # Reset boundary integrations
    intFaceLD2 = 0.0
    intFaceRD2 = 0.0

    # Loop over points of Domain 2
    c_0 = 0
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                # Loop over left boundary points of Domain 2
                if ((points[c_0][0]-x0) <= tolerance):
                    points_fetch = [points[c_0][1], points[c_0][2]]
                    # Fetch data from the other solver
                    scalar_field_A[c_0] = MUI_Interfaces["interface2D01"].\
                                                        fetch(name_fetchA,
                                                        points_fetch,
                                                        t,
                                                        s_sampler,
                                                        t_sampler)

                    # Calculate the integration of left boundary points of Domain 2
                    intFaceLD2 += scalar_field_A[c_0]
                else: # Loop over 'internal' points of Domain 2

                    #Calculate the diffusion of pseudo scalar field of Domain 2
                    scalar_field_A[c_0] += dr*(scalar_field_A[c_0 - 1] - scalar_field_A[c_0])

                    # Loop over right boundary points of Domain 2
                    if (( (x0+lx) - points[c_0][0]) <= (tolerance)):
                        points_push = [(y0+(ly/(Ny-1))*j), (z0pr+(lzpr/(Nz-1))*k)]
                        # push data to the other solver
                        MUI_Interfaces["interface2D02"].push(name_pushA, points_push, scalar_field_A[c_0])
                        # Calculate the integration of right boundary points of Domain 2
                        intFaceRD2 += scalar_field_A[c_0]
                c_0 += 1

    # Commit 't' step of MUI
    MUI_Interfaces["interface2D02"].commit(t)

    # Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
    if ( (t - forgetSteps) > 0 ):
        MUI_Interfaces["interface2D02"].forget(t - forgetSteps)

    # Output the pseudo scalar field and the boundary integrations
    if ((t % outputInterval) == 0):
        outputFileName = "coupling_results" + str(rank) + "/scalar_field_coarse_" + str(t) + ".csv"
        outputFile = open(outputFileName,"w+")
        outputFile.write("\"X\",\"Y\",\"Z\",\"scalar_field_A\"\n")
        c_0 = 0
        for k in range(Nz):
            for j in range(Ny):
                for i in range(Nx):
                    outputFile.write("%f,%f,%f,%f\n" % (points[c_0][0],points[c_0][1],points[c_0][2],scalar_field_A[c_0]))
                    c_0 += 1

        outputFile.close() 

        outputIntegration = open(outputIntegrationName,"a")
        outputIntegration.write("%i,%f,%f\n" % (t,intFaceLD2,intFaceRD2))
        outputIntegration.close()
