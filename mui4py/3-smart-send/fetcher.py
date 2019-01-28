import mui4py
import time

appcomm = mui4py.mpi_split_by_app()
rank = appcomm.Get_rank()
domain = "fetcher"

# Configuration for the interfaces
config = mui4py.Config(2, mui4py.FLOAT64)
mui4py.set_default_config(config)
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI)
uniface.set_data_types({"data": mui4py.STRING})
t_sampler = mui4py.ChronoSamplerExact
s_sampler = mui4py.SamplerExact

# Rank to geometry mapping
if rank == 0:
    box_p1 = [-1.0, -1.0]
    box_p2 = [0.0, 0.0]
elif rank == 1:
    box_p1 = [0.0, -1.0]
    box_p2 = [1.0, 0.0]
elif rank == 2:
    box_p1 = [-1.0, 0.0]
    box_p2 = [0.0, 1.0]
elif rank == 3:
    box_p1 = [0.0, 0.0]
    box_p2 = [1.0, 1.0]

# Point at the center of the quadrant
fetch_point = [(box_p2[0]+box_p1[0])/2.0, (box_p2[1]+box_p1[1])/2.0]

# Initialisation of smart sending
recv_box = mui4py.geometry.Box(box_p1, box_p2)
uniface.announce_recv_span(0, 4, recv_box) 

steps = 4
sleep_t = 0.1
appcomm.Barrier()
time.sleep(sleep_t)
for s in range(1,steps+1):
    if rank == 0:
        print("\nStep ", s, flush=True)
    # Let time to print
    time.sleep(sleep_t)
    output = ""
    fetch_val = uniface.fetch("data", fetch_point, s, s_sampler(), t_sampler())
    output += "rank: " +  str(rank) + " fetching at " + str(fetch_point) + " (step: %d)" % s
    print (output, "\n", "   Domain:", domain, "Fetched value: '{}'".format(fetch_val), flush=True)
    # Let time to print
    time.sleep(sleep_t)
