import mui4py
import mpi4py
import sys
import time
import numpy as np

domain = "pusher"
push_val = 666.0

# Get intra communicator
appcomm = mui4py.mpi_split_by_app()

# Do not print Warnings or Info messages to not interfere
# with profiling.
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

# Iterations
Nt = int(sys.argv[1])
Ni = int(sys.argv[2])
Nj = int(sys.argv[3])
Npoints = Ni*Nj

# Allocate memory for the points. Speedup things.
uniface.push_size("data_py", Npoints, False)
uniface.push_size("data_many", Npoints, False)
uniface.push_size("data_raw", Npoints, False)

# Setup data. Random array from predefined seed.
np.random.seed(666)
vals = np.random.rand(Npoints)
points = np.zeros((Npoints, 2))
c = 0
# Create array of points
for i in range(Ni):
    for j in range(Nj):
        points[c] = [i, j]
        c += 1

## USING PYTHONIC INTERFACE

t_push_py = 0.0
t_push_py_max = 0.0
t_push_py_min = 100.0

# Time loop
for t in range(0, Nt):
    for i, p in enumerate(points):
        t0 = time.time()
        uniface.push("data_py", p, vals[i])
        t_push_py += time.time() - t0
        t_new = time.time() - t0
        t_push_py_max = max(t_new, t_push_py_max)
        t_push_py_min = min(t_new, t_push_py_min)
        t_push_py += t_new

    uniface.commit(t)

## USING PUSH_MANY

t_push_many = 0.0

# Time loop
for t in range(Nt,  2*Nt):
    t0 = time.time()
    uniface.push_many("data_many", points, vals)
    t_push_many += time.time() - t0
    uniface.commit(t)


## USING RAW INTERFACE

t_push_raw = 0.0
t_push_raw_max = 0.0
t_push_raw_min = 100.0

# Time loop
for t in range(2*Nt, 3*Nt):
    for i, p in enumerate(points):
        t0 = time.time()
        uniface.raw.push_double("data_raw", uniface.raw_point(p), vals[i])
        t_push_raw += time.time() - t0
        t_new = time.time() - t0
        t_push_raw_max = max(t_new, t_push_raw_max)
        t_push_raw_min = min(t_new, t_push_raw_min)
        t_push_raw += t_new
    uniface.commit(t)

## Print timing results

mpi4py.MPI.COMM_WORLD.Barrier()
time.sleep(1)
print("\nPush function times in seconds:", flush=True)
nof_push_calls = Nt * Ni * Nj 
print("\tPython total time: {:.10f}\
       \n\tPython min time: {:.2E}\
      \n\tPython max time: {:.2E}\
      \n\tPython mean time: {:.2E}\
      \n\tRaw total time: {:.10f}\
      \n\tRaw min time: {:.2E}\
      \n\tRaw max time: {:.2E}\
      \n\tRaw mean time: {:.2E}\
      \n\tPush_many total time: {:.10f}\
      \n--> Raw 'push' call was {:.10f} times faster than Pythonic.\
      \n--> 'push_many' call was {:.10f} times faster than Pythonic."\
                                    .format(t_push_py, t_push_py_min, t_push_py_max,
                                            t_push_py/nof_push_calls, t_push_raw,
                                            t_push_raw_min, t_push_raw_max,
                                            t_push_raw/nof_push_calls, t_push_many,
                                            t_push_py/t_push_raw, t_push_raw/t_push_many), flush=True)
