import sys
import mui4py
domain = sys.argv[1]
push_val = 666.0
push_val_assign = 333.0

# Get intra communicator
appcomm = mui4py.mpi_split_by_app()
rank = appcomm.Get_rank()

print("\nInfo:")
print("\tMpi version:", mui4py.get_mpi_version())
print("\tCompiler version:", mui4py.get_compiler_version())
print("\tCompiler config:", mui4py.get_compiler_config())

# Configuration for the interface (dimensions, float type, int type)
dims = 2
config = mui4py.Config(dims, mui4py.FLOAT64)

# Default configuration for every object so no need to pass it as an argument 
# to any class constructor.
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI, config=config)
uniface.set_data_types({"data": mui4py.FLOAT,
                       "data_assign": mui4py.FLOAT})

# Push/fetch some data
point1 = [1.0, 1.0]
point2 = [2.0, 1.0]
t = 0.0
# Value at a certain point
uniface.push("data", point1, push_val)
uniface.push("data", point2, push_val)

# Value spatially independent
uniface.assign("data_assign", push_val_assign)
uniface.commit(t)
