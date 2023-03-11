import numpy as np
import sys
import time
import mui4py
from mpi4py import MPI
domain = sys.argv[1]

# Configuration for the interfaces
config1d = mui4py.Config(1, mui4py.FLOAT64)
config2d = mui4py.Config(2, mui4py.FLOAT32)
config3d = mui4py.Config(3, mui4py.FLOAT64)

appcomm = mui4py.mpi_split_by_app()
rank = appcomm.Get_rank()

# Set up interfaces by type
# 1D interfaces
ifs = []
ifaces1d = mui4py.create_unifaces(domain, ifs, config1d, MPI.COMM_WORLD)

# 2D interfaces
ifs = ["I1", "I2"]
ifaces2d = mui4py.create_unifaces(domain, ifs, config2d) 
ifaces2d["I1"].set_data_types({"data1": mui4py.FLOAT64,
                              "data2": mui4py.STRING})
ifaces2d["I2"].set_data_types({"data": mui4py.INT64})

# 3D interfaces
ifs = ["I3"]
ifaces3d = mui4py.create_unifaces(domain, ifs, config3d) 
ifaces3d["I3"].set_data_types({"more_data": mui4py.INT64})

time.sleep(0.1)
print ("Domain", domain + ",", "rank", rank, "connected to interfaces:", 
        "\n\t1D:", list(ifaces1d.keys()),
        "\n\t2D:", list(ifaces2d.keys()), 
        "\n\t3D:", list(ifaces3d.keys()))
