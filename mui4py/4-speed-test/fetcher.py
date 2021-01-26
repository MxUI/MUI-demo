import mpi4py
import mui4py
import sys
import time
import numpy as np

domain = "fetcher"

# Get intra communicator
appcomm = mui4py.mpi_split_by_app()

# Do not print Warnings or Info messages to not interfere
# with profiling
mui4py.set_quiet(True)

# Configuration for the interface (dimensions, float type, int type)
dims = 2
config = mui4py.Config(dims, mui4py.FLOAT64)

# Default configuration for every object so no need to pass it as an argument 
# to any class constructor.
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI, config=config)
uniface.set_data_types({"data_py": mui4py.FLOAT64,
                        "data_many": mui4py.FLOAT64,
                        "data_raw": mui4py.FLOAT64})

# Push/fetch point

# Iterations
Nt = int(sys.argv[1])
Ni = int(sys.argv[2])
Nj = int(sys.argv[3])
Npoints = Ni*Nj

# Setup data. Random array from predefined seed.
np.random.seed(666)
push_vals = np.random.rand(Npoints)
points = np.zeros((Npoints, 2))
c = 0
for i in range(Ni):
    for j in range(Nj):
        points[c] = [i, j]
        c += 1

print("Running pythonic...", flush=True)

## -- USING PYTHONIC INTERFACE --

# Spatial and temporal samplers
s_sampler_exact = mui4py.SamplerExact()
t_sampler_exact = mui4py.ChronoSamplerExact()

t_fetch_py = 0.0
t_fetch_py_max = 0.0
t_fetch_py_min = 100.0

fetch_vals = np.zeros(Npoints)
# Time loop
for t in range(0, Nt):
    # For i, j point in the grid
    for i, p in enumerate(points):
        t0 = time.time()
        fetch_vals[i]  = uniface.fetch("data_py", p, t,
                                       s_sampler_exact, t_sampler_exact)
        t_fetch_py += time.time() - t0
        t_new = time.time() - t0
        t_fetch_py_max = max(t_new, t_fetch_py_max)
        t_fetch_py_min = min(t_new, t_fetch_py_min)
        t_fetch_py += t_new

if (fetch_vals == push_vals).all():
    print("Sucess fetch_py!", flush=True)

print("Running fetch/push_many...",flush=True)
t_fetch_many = 0.0

fetch_vals = None
# # Time loop
for t in range(Nt, 2*Nt):
    t0 = time.time()
    fetch_vals = uniface.fetch_many("data_many", points, t,
                                    s_sampler_exact, t_sampler_exact)
    t_fetch_many += time.time() - t0

if (fetch_vals == push_vals).all():
    print("Sucess fetch_many!", flush=True)

## USING RAW INTERFACE
print("Running raw...",flush=True)

t_fetch_raw = 0.0
t_fetch_raw_max = 0.0
t_fetch_raw_min = 100.0

# Spatial and temporal samplers manually configured
# NOTE: Using the raw interface one has to manually
#       define and configure one sampler per data type and 
#       specific configuration of the interface.
s_sampler_raw = mui4py.SamplerExact()
s_sampler_raw.configure(config, io_data_type=mui4py.FLOAT64)
t_sampler_raw = mui4py.ChronoSamplerExact()
t_sampler_raw.configure(config)

fetch_vals = np.zeros(Npoints)
# Time loop
for t in range(2*Nt, 3*Nt):
    for i, p in enumerate(points):
        t0 = time.time()
        fetch = uniface.raw.fetch_double_exact_chrono_exact
        fetch_vals[i] = fetch("data_raw", uniface.raw_point(p), t,
                          s_sampler_raw.raw, t_sampler_raw.raw, True)
        t_new = time.time() - t0
        t_fetch_raw_max = max(t_new, t_fetch_raw_max)
        t_fetch_raw_min = min(t_new, t_fetch_raw_min)
        t_fetch_raw += t_new

if (fetch_vals == push_vals).all():
    print("Sucess fetch_raw!", flush=True)

mpi4py.MPI.COMM_WORLD.Barrier()
time.sleep(1)
print("\nFetch function times in seconds:", flush=True)
nof_fetch_calls = Nt * Ni * Nj
print("\tPython total time: {:.10f}\
       \n\tPython min time: {:.2E}\
       \n\tPython max time: {:.2E}\
       \n\tPython mean time: {:.2E}\
       \n\tRaw total time: {:.10f}\
       \n\tRaw min time: {:.2E}\
       \n\tRaw max time: {:.2E}\
       \n\tRaw mean time: {:.2E}\
       \n\tFetch_many total time: {:.10f}\
       \n-> Raw 'fetch' call was {:.10f} times faster than Pythonic.\
       \n-> 'fetch_many' call was {:.10f} times faster than Raw."\
                            .format(t_fetch_py, t_fetch_py_min, t_fetch_py_max,
                                    t_fetch_py/nof_fetch_calls, t_fetch_raw,
                                    t_fetch_raw_min, t_fetch_raw_max,
                                    t_fetch_raw/nof_fetch_calls, t_fetch_many,
                                    t_fetch_py/t_fetch_raw, t_fetch_raw/t_fetch_many), flush=True)
