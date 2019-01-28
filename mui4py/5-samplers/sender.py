import mui4py
from mpi4py import MPI
import sys
import time
domain = sys.argv[1]

appcomm = mui4py.mpi_split_by_app()

# Configuration for the interface
dims = 3
config = mui4py.Config(dims, mui4py.FLOAT64, mui4py.INT32, force_casting=True)
mui4py.set_default_config(config)

# Interface 
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI)

# Set types for data tags in push/fetch calls
uniface.set_data_types({ "data": mui4py.FLOAT64})

# Push/fetch points to a 4x4 grid
# 1 |  1 2 3 4
# 2 |  4 2 3 1
# 3 |  1 0 1 0
# 4 |  0 1 0 1
#    ---------- 
#      1 2 3 4

point = [0.5, 0.5, 0.5]

# EXAMPLE 1
uniface.push("data", point, 
             "A boring string from domain: %s" % domain)
uniface.commit(0)
fetch_val = uniface.fetch("data_string", point, 0, s_sampler_exact,
                          t_sampler_exact)

# EXAMPLE 2: 
