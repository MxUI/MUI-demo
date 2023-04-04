"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library Demo 10-2-3           #
#                                                                            #
# Copyright (C) 2023 E. R. Fernandez                                         #
#                                                                            #
# This software is jointly licensed under the Apache License, Version 2.0    #
# and the GNU General Public License version 3, you may use it according     #
# to either.                                                                 #
#                                                                            #
# ** Apache License, version 2.0 **                                          #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# You may obtain a copy of the License at                                    #
#                                                                            #
# http://www.apache.org/licenses/LICENSE-2.0                                 #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
# ** GNU General Public License, version 3 **                                #
#                                                                            #
# This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by       #
# the Free Software Foundation, either version 3 of the License, or          #
# (at your option) any later version.                                        #
#                                                                            #
# This program is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
# GNU General Public License for more details.                               #
#                                                                            #
# You should have received a copy of the GNU General Public License          #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
##############################################################################
#
# @file fetcher.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Fetcher file of MUI Python wrapper smart send demo.
#
"""

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
t_sampler = mui4py.TemporalSamplerExact
s_sampler = mui4py.SamplerExact
synchronised=False
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
uniface.announce_recv_span(0, 4, recv_box, synchronised)

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
