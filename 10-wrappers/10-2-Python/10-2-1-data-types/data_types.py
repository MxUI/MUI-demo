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
# @file data_types.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Python file of MUI Python wrapper data types demo.
#
"""

import mui4py
from mpi4py import MPI
import sys
import numpy as np
import time
domain = sys.argv[1]

appcomm = mui4py.mpi_split_by_app()

# Configuration for the interface
dims = 3
# Option force_casting=True, forces the casting of data pushed/fetched
# to the type specified with set_data_types(). This cannot be changed later after
# an object has been created with such configuration.
config = mui4py.Config(dims, mui4py.FLOAT64)
# Default configuration for every object. No need to pass it as an argument 
# to any class constructor.
mui4py.set_default_config(config)

# Interface 
URI = "mpi://" + domain + "/ifs1"
uniface = mui4py.Uniface(uri=URI)

#NOTE 1: 'mui4py.FLOAT' will match native python 'int' type and
#        convert it appropriately to either 'mui4py.FLOAT64' 
#        or 'mui4py.FLOAT32'. Same applies for 'mui4py.INT'.
#NOTE 2: mui4py.[INT|FLOAT][32|64] are aliases for numpy types.

# Set types for data tags in push/fetch calls
uniface.set_data_types({"data_string": mui4py.STRING,
                        "data_float": mui4py.FLOAT,
                        "data_float64": mui4py.FLOAT64,
                        "data_float32": mui4py.FLOAT32,
                        "data_int": mui4py.INT,
                        "data_int32": mui4py.INT32,
                        "data_int64": mui4py.INT64})

# Spatial and time samplers. The appropriate type dispatching for the
# samplers is done withing fetch/push using 'Config()' and 'set_data_type()' info.
t_sampler_exact = mui4py.TemporalSamplerExact()
s_sampler_exact = mui4py.SamplerExact()

# Header/footer printing syncronized functions for the examples
def print_header(dom, example):
    MPI.COMM_WORLD.Barrier()
    if dom == "dom1":
        print("\nEXAMPLE %d." % example, flush=True)
        time.sleep(0.1)
    MPI.COMM_WORLD.Barrier()

def print_footer(push_val, fetch_val, domain):
    MPI.COMM_WORLD.Barrier()
    print("\tPushed '%s' in domain %s." % (push_val, domain), flush=True)
    print("\tFetched '%s' in domain %s." % (fetch_val, domain), flush=True)
    time.sleep(0.1)
    MPI.COMM_WORLD.Barrier()

# Push/fetch point
point = [0.5, 0.5, 0.5]

# EXAMPLE 1
print_header(domain, 1)
t = 0
push_val = "A boring string from domain: %s" % domain
uniface.push("data_string", point, push_val)
uniface.commit(t)
fetch_val = uniface.fetch("data_string", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_val, fetch_val, domain)

# EXAMPLE 2: 
#   - Forced conversion from int value = 10 to float32
#   - No harm as that value is withing a range it can be safely represented by a float32
print_header(domain, 2)
t = 1
push_val = 10
uniface.push("data_float32", point, push_val)
uniface.commit(t)
fetch_val = uniface.fetch("data_float32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_val, fetch_val, domain)

# EXAMPLE 3: 
#   - Forced conversion from int value = 100000 to float32
#   - This is not allowed since the value cannot be safely represented.
print_header(domain, 3)
t = 2
push_val = 100000
try:
    uniface.push("data_float32", point, push_val)
except Exception as err:
    print("Error: ", err)
uniface.commit(t)
fetch_val = uniface.fetch("data_float32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_val, fetch_val, domain)

# EXAMPLE 4: 
#   - Forced conversion from string value = "23" to int
#   - Conversion from string to numeric types allowed.
print_header(domain, 4)
t = 3
push_val = "23"
uniface.push("data_int32", point, push_val)
uniface.commit(t)
fetch_val = uniface.fetch("data_int32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_val, fetch_val, domain)

# EXAMPLE 5: 
#   - Forced conversion from int64 to int32 with a value that causes overflow.
#   - This is not allowed
print_header(domain, 5)
t = 4
push_value = np.int64(100000000000)
try:
    uniface.push("data_int32", point, push_value)
except Exception as err:
    print("Error: ", err)
uniface.commit(t)
fetch_val = uniface.fetch("data_int32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_value, fetch_val, domain)

# EXAMPLE 6: 
#   - Forced conversion from int64 to int32 with a value that do NOT overflow.
#   - This is allowed.
print_header(domain, 6)
t = 5
push_val  = np.int64(1000)
try:
    uniface.push("data_int32", point, push_val)
except Exception as err:
    print("Error: ", err)
uniface.commit(t)
fetch_val = uniface.fetch("data_int32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_val, fetch_val, domain)

# EXAMPLE 7: 
#   - Forced conversion from float value = 10.2 to int
#   - Truncation happens so an exception is thrown as it is not safe.
print_header(domain, 7)
t = 6
push_value = 10.2
try:
    uniface.push("data_int32", point, push_value)
except Exception as err:
    print("Error: ", err)
uniface.commit(t)
fetch_val = uniface.fetch("data_int32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_value, fetch_val, domain)

# EXAMPLE 8:
#   - Forced conversion from float64 value = 1.1231321355 to float32
#   - The decimal part is truncated to float32 precision. This is legal but use with caution! 
print_header(domain, 8)
t = 7
push_val = np.float64(1.1231321355)
uniface.push("data_float32", point, push_val)
uniface.commit(t)
fetch_val = uniface.fetch("data_float32", point, t,
                          s_sampler_exact,
                          t_sampler_exact)
print_footer(push_val, fetch_val, domain)
