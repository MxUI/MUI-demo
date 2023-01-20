"""
#
#% @file Pusher_Fetcher_1.py
#% @author W.L.
#% @date 11 Jan 2019
#% @brief Pusher_Fetcher_1.py <-> MUI (Python) <-> MUI (C++) <-> Pusher_Fetcher_0.cpp two way Coupling test case.
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
			points_push[c_0] = [i+4.0, j*0.1, k+0.0]
			c_0 += 1

push_values = np.zeros(Npoints)
for i in range(Npoints):
	push_values[i] = 48.0

# define fetch points and evaluation
points_fetch = np.zeros((Npoints, dimensionMUI))
c_1 = 0
for i in range(Ni):
	for j in range(Nj):
		for k in range(Nk):
			points_fetch[c_1] = [i+6.5, j*0.1, k+0.0]
			c_1 += 1		

fetch_vals = np.zeros(Npoints)

# Spatial/temporal samplers
t_sampler = mui4py.ChronoSamplerExact()
s_sampler = mui4py.SamplerPseudoNearestNeighbor(r)

# # Fetch ZERO step
iface.barrier(0)

for n in range(1, steps):

	if rank == 0:
		print("\n{PUSHER_FETCHER_1} Step ", n, flush=True)

	# MUI Push boundary points and commit current steps
	for i, p in enumerate(points_push):
		iface.push("data_python", p, push_values[i])
	iface.commit(n)

    # MUI Fetch internal points
	fetch_vals  = iface.fetch_many("data_cpp", points_fetch, n,
									s_sampler, 
									t_sampler)

	if (rank == 0):
		print ("{PUSHER_FETCHER_1} cpp_fetch: ", fetch_vals)