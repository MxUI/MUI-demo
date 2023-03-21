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
# @file heat-right.py
# @author W. Liu
# @date 19 March 2019
# @brief Right domain on Aitken's dynamic points to demonstrate MUI coupling algorithm.
#
"""

from mpi4py import MPI
import mui4py
import numpy as np
import csv
import os

#                             Left Domain                       Right Domain             
# Coarse : +-------+-------+-------+-------o=======+=======o-------+-------+-------+-------+
#          0       1       2       3       4       5       6       7       8       9      10
# +: grid points
# o: interface points
# -: single domain zone
# =: overlapping zone

# MUI parameters
dimensionMUI = 1
data_types = {"u0": mui4py.FLOAT64,
                "u": mui4py.FLOAT64}

# MUI interface creation
config3d = mui4py.Config(dimensionMUI, mui4py.FLOAT64)
URI = "mpi://right/ifs" 
iface = mui4py.Uniface(uri=URI, config=config3d) 
iface.set_data_types(data_types)

# App common world claims 
MPI_COMM_WORLD = mui4py.mpi_split_by_app()

# Local communicator rank and size
rank = MPI_COMM_WORLD.Get_rank()
size = MPI_COMM_WORLD.Get_size()

u1 = np.zeros(110)
u2 = np.zeros(110)

# Initialise values from file
content = []
try:
    with open('Resources/right_AITKEN.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            content.append(row)
        for i in range(40, 110, 10):
            u1[i] = float(content[int(i*0.1-3)][1])
except FileNotFoundError:
    print("right_AITKEN.csv missing")
    for i in range(40, 110, 10):
        u1[i] = float(0.0)

# Create folder if don't exists
makedirMString = "results_right" + str(rank)
if not os.path.exists(makedirMString):
    os.makedirs(makedirMString)
    print("Folder " , makedirMString, " created ")
else:
    print("Folder " , makedirMString,  " already exists")
makedirMIterString = "results_iteration_right" + str(rank)
if not os.path.exists(makedirMIterString):
    os.makedirs(makedirMIterString)
    print("Folder " , makedirMIterString, " created ")
else:
    print("Folder " , makedirMIterString,  " already exists")

k, H = 0.515, 1.0

ptsVluInit = []

for i in range(40, 110, 10):
    pt = iface.Point([i])
    ptsVluInit.append((pt, u1[i]))

# Samplers
t_sampler = mui4py.TemporalSamplerExact()
s_sampler = mui4py.SamplerPseudoNearestNeighbor(30)
a_algorithm = mui4py.AlgorithmAitken(0.01, 1.0, MPI_COMM_WORLD, ptsVluInit, 0.0, config=config3d)

# Print off a hello world message
print("Hello world from Right rank ", rank, " out of ", size, " MUI processors\n")

# Output
filenameR = f"results_right{rank}/solution-right_AITKEN_0.csv"
with open(filenameR, "w") as outputFileRight:
    outputFileRight.write("\"X\",\"u\"\n")
    for i in range(40, 110, 10):
        outputFileRight.write(f"{i*H},{u1[i]},\n")
filenameIterR = f"results_iteration_right{rank}/solution-right_AITKEN_0.csv"
with open(filenameIterR, "w") as outputFileIterRight:
    outputFileIterRight.write("\"X\",\"u\"\n")
    for i in range(40, 110, 10):
        outputFileIterRight.write(f"{i*H},{u1[i]},\n")

for t in range(1, 11):
    for iter in range(1, 101):
        print(f"Right grid time {t} iteration {iter}")
        u1[40] = iface.fetch("u", 40 * H, t, iter, s_sampler, t_sampler, a_algorithm)
        
        if ((t>=4) and (t<6)):
            u1[42] = iface.fetch("u", 42 * H, t, iter, s_sampler, t_sampler, a_algorithm)

        print(f"Right under relaxation factor at t= {t} iter= {iter} is {a_algorithm.get_under_relaxation_factor(t, iter)}")
        print(f"Right residual L2 Norm at t= {t} iter= {iter} is {a_algorithm.get_residual_L2_Norm(t, iter)}")

        # calculate 'interior' points
        for i in range(50, 100, 10):
            u2[i] = u1[i] + k / (H * H) * (u1[i - 10] + u1[i + 10] - 2 * u1[i])
        # calculate 'boundary' points
        u2[100] = 0.0
        u2[40] = u1[40]

        if ((t>=4) and (t<6)):
            u2[42] = u1[42]

        # push data to the other solver
        iface.push("u0", 60 * H, u1[60])
        iface.commit(t, iter)

        # I/O
        u1, u2 = u2, u1

        # Output
        filenameR = f"results_iteration_right{rank}/solution-right_AITKEN_{((t-1)*100) + iter}.csv"
        with open(filenameR, "w") as outputFileRight:
            outputFileRight.write("\"X\",\"u\"\n")
            outputFileRight.write(f"{40*H},{u1[40]},\n")
            if ((t>=4) and (t<6)):
                outputFileRight.write(f"{42*H},{u1[42]},\n")
            for i in range(50, 110, 10):
                outputFileRight.write(f"{i*H},{u1[i]},\n")

    # Output
    filenameR = f"results_right{rank}/solution-right_AITKEN_{t}.csv"
    with open(filenameR, "w") as outputFileRight:
        outputFileRight.write("\"X\",\"u\"\n")
        outputFileRight.write(f"{40*H},{u1[40]},\n")
        if ((t>=4) and (t<6)):
            outputFileRight.write(f"{42*H},{u1[42]},\n")
        for i in range(50, 110, 10):
            outputFileRight.write(f"{i*H},{u1[i]},\n")
