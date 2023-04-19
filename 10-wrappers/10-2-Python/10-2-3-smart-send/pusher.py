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
# @file pusher.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Pusher file of MUI Python wrapper smart send demo.
#
"""

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
