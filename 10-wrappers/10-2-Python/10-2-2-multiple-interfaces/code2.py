import numpy as np
import time
import sys
import mui4py
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
ifaces1d = mui4py.create_unifaces(domain, ifs, config1d) 

# 2D interfaces
ifs = []
ifaces2d = mui4py.create_unifaces(domain, ifs, config2d) 
# 3D interfaces

ifs = ["I3"]
ifaces3d = mui4py.create_unifaces(domain, ifs, config3d) 
ifaces3d["I3"].set_data_types({"more_data": mui4py.INT64})

# Join all interfaces under a single dict for convenience
ifaces = {**ifaces1d, **ifaces2d, **ifaces3d}
time.sleep(0.1)
print ("Domain", domain + ",", "rank", rank, "connected to interfaces:", 
        "\n\t1D:", list(ifaces1d.keys()),
        "\n\t2D:", list(ifaces2d.keys()), 
        "\n\t3D:", list(ifaces3d.keys()))
