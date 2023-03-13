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
# @file pusher.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Pusher of MUI Python wrapper simple demo.
#
"""

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
uniface.push("data_assign", push_val_assign)
uniface.commit(t)
