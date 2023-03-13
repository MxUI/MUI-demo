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
# @file fetcher.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Fetcher of MUI Python wrapper simple demo.
#
"""

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
