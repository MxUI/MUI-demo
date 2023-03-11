import sys
import mui4py
domain = sys.argv[1]

# Get intra communicator
appcomm = mui4py.mpi_split_by_app()
rank = appcomm.Get_rank()

# Configuration for the interface (dimensions, float type, int type)
dims = 2
config = mui4py.Config(dims, mui4py.FLOAT64)

# Default configuration for every object so no need to pass it as an argument 
# to any class constructor.
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI, config=config)
uniface.set_data_types({"data": mui4py.FLOAT,
                       "data_assign": mui4py.FLOAT})

# Spatial and temporal samplers
t_sampler_exact = mui4py.TemporalSamplerExact()
s_sampler_exact = mui4py.SamplerExact()

# Push/fetch some data
point = [1.0, 1.0]
t = 0.0

# Wait until commit has completed in 'pusher'
uniface.barrier(0)
fetch_val_assign = uniface.fetch("data_assign")
print("\nFetched spatially independent value '{}' from domain {}.".format(fetch_val_assign, domain))
fetch_val = uniface.fetch("data", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
#print ("Points available:", uniface.fetch_points("data", t))
print("Fetched spatially dependent value '{}' from domain {}.".format(fetch_val, domain))
