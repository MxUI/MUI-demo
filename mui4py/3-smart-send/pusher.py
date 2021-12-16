import mui4py
domain = "pusher"

# Decomposition of the domain [-1,-1] x [1,1] for receiver  processes
# x: Fetched points at the center of quadrants
# Ri: Ranks
#
#  (-1,1) |=====|=====| (1,1)
#         | R2  | R3  |
#         |  x  |  x  |
#         |     |     |
#         |=====|=====|
#         | R0  | R1  |
#         |  x  |  x  |
#         |     |     |
# (-1,-1) |=====|=====| (1,-1)
#    

appcomm = mui4py.mpi_split_by_app()

# Configuration for the interfaces
config = mui4py.Config(2, mui4py.FLOAT64)
# Default configuration for every object.
mui4py.set_default_config(config)
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI)
uniface.set_data_types({"data": mui4py.STRING})
synchronised=False
# Eeach step from 1 to 3 sender uses a different span:
#     s = 1 -> Box intersecting R0 and R2
#     s = 2 -> Sphere intersecting R1 and R3
#     s = 3 -> Point intersecting R3
#     s = 4 -> Smart sending is disabled due to timeout
send_region_t1 = mui4py.geometry.Box([-1, -1], [-0.5, 1.0])
send_region_t2 = mui4py.geometry.Sphere([1.0, 0.0], 0.5)
send_region_t3 = mui4py.geometry.Point([0.5, 0.5])

# Initialisation of smart sending
uniface.announce_send_span(0, 1, send_region_t1,synchronised)

steps = 4
# At s=0  (Smart sending do not work in the first step)
print ("\nLooping...", flush=True)
for s in range(1,steps+1):
    # Change send span
    if s == 1:
        uniface.announce_send_span(2, 2, send_region_t2,synchronised)
    elif s == 2:
        uniface.announce_send_span(3, 3, send_region_t3,synchronised)
    uniface.push("data", [-0.5, -0.5], "R0_t%d" % s)
    uniface.push("data", [0.5, -0.5], "R1_t%d" % s)
    uniface.push("data", [0.5, 0.5], "R2_t%d" % s)
    uniface.push("data", [-0.5, 0.5], "R3_t%d" % s)
    uniface.commit(s)
