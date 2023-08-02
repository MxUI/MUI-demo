"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library Demo 10-2-5           #
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
# @file Pusher_Fetcher_1.py
# @author W. Liu
# @date 11 January 2019
# @brief Pusher_Fetcher_1.py <-> MUI (Python) <-> MUI (C++) <-> Pusher_Fetcher_0.cpp two way Coupling test case.
#
"""

from mpi4py import MPI
import datetime
import mui4py
import numpy as np
import time

# MUI parameters
dimensionMUI = 3
data_types = {"data_python": mui4py.FLOAT64,
				"data_cpp": mui4py.FLOAT64}

# MUI interface creation
config3d = mui4py.Config(dimensionMUI, mui4py.FLOAT64)
URI = "mpi://PUSHER_FETCHER_1/COUPLING" 
iface = mui4py.Uniface(uri=URI, config=config3d) 
iface.set_data_types(data_types)

# App common world claims 
MPI_COMM_WORLD = mui4py.mpi_split_by_app()

# Local communicator rank
rank = MPI_COMM_WORLD.Get_rank()

steps = 10				# number of time steps
iterations = 10			# number of iterations per step
r = 1.					# search radius
Ni = int(1)
Nj = int(11)
Nk = int(1)
Npoints = Ni*Nj*Nk

# define push points and evaluation
points_push = np.zeros((Npoints, dimensionMUI))
c_0 = 0
for i in range(Ni):
	for j in range(Nj):
		for k in range(Nk):
			points_push[c_0] = [i+4.0+rank, j*0.1, k+0.0]
			c_0 += 1

# define fetch points
fetch_points_only = np.zeros((Npoints, dimensionMUI))

push_values = np.zeros(Npoints)
for i in range(Npoints):
	push_values[i] = 48.0

# define fetch points and evaluation
points_fetch = np.zeros((Npoints, dimensionMUI))
c_1 = 0
for i in range(Ni):
	for j in range(Nj):
		for k in range(Nk):
			points_fetch[c_1] = [i+6.0+rank, j*0.1, k+0.0]
			c_1 += 1		

fetch_vals = np.zeros(Npoints)
fetch_vals_only = np.zeros(Npoints)

# Define and announce MUI send/receive span
send_span = mui4py.geometry.Box([i+4.0+rank, 0, 0], [i+4.0+rank, 20, 1])
recv_span = mui4py.geometry.Box([i+6.0+rank, 0, 0], [i+6.0+rank, 20, 1])
iface.announce_recv_span(1, (steps), recv_span, False)
iface.announce_send_span(1, (steps), send_span, False)

# Spatial/temporal samplers
t_sampler = mui4py.TemporalSamplerExact()
s_sampler = mui4py.SamplerPseudoNearestNeighbor(r)

# # Fetch ZERO step
iface.barrier(0, 0)

for n in range(1, steps):
	for iter in range(1, iterations):
		if rank == 0:
			print("\n{PUSHER_FETCHER_1} Step ", n, "Iteration ", iter, flush=True)

		# MUI Push boundary points and commit current steps -- both `push()` and `push_many()` have been tested here
		if (n < 5):
			for i, p in enumerate(points_push):
				iface.push("data_python", p, push_values[i])
		else:
			iface.push_many("data_python", points_push, push_values)

		commit_return = iface.commit(n,iter)
		if (rank == 0):
			print ("{PUSHER_FETCHER_1} commit_return: ", commit_return)

		# MUI is_ready function
		is_ready_return = iface.is_ready("data_cpp", n, iter)
		if (rank == 0):
			print ("{PUSHER_FETCHER_1} is_ready_return: ", is_ready_return)

		# MUI Fetch internal points -- both `fetch()` and `fetch_many()` have been tested here
		if (n < 5):
			for i, p in enumerate(fetch_vals):
				fetch_vals[i]  = iface.fetch("data_cpp", points_fetch[i], n, iter,
												s_sampler,
												t_sampler)
		else:
			fetch_vals  = iface.fetch_many("data_cpp", points_fetch, n, iter,
											s_sampler,
											t_sampler)

		# MUI forget function
		if (n > 2):
			t = (n-2, iter)
			iface.forget(t)

		# MUI set_memory function
		iface.set_memory(2)

		# Print fetched values
		if (rank == 0):
			print ("{PUSHER_FETCHER_1} cpp_fetch: ", fetch_vals)