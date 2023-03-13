"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library                       #
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
# @file code3.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief code3 file of MUI Python wrapper multiple interfaces demo.
#
"""

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
ifs = ["I4"]
ifaces1d = mui4py.create_unifaces(domain, ifs, config1d) 
ifaces1d["I4"].set_data_types({"data": mui4py.INT32})
# 2D interfaces
ifs = []
ifaces2d = mui4py.create_unifaces(domain, ifs, config2d) 
# 3D interfaces
ifs = []
ifaces3d = mui4py.create_unifaces(domain, ifs, config3d) 

# Join all interfaces under a single dict for convenience
ifaces = {**ifaces1d, **ifaces2d, **ifaces3d}
time.sleep(0.1)
print ("Domain", domain + ",", "rank", rank, "connected to interfaces:", 
        "\n\t1D:", list(ifaces1d.keys()),
        "\n\t2D:", list(ifaces2d.keys()), 
        "\n\t3D:", list(ifaces3d.keys()))
