/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 8                *
*                                                                            *
* Copyright (C) 2020 S. M. Longshaw                                          *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
******************************************************************************/

/**
 * @file fetchall.cpp
 * @author S. M. Longshaw
 * @date 25 September 2020
 * @brief Simple "Hello World" example demonstrating the non-interpolating
 *        fetch_points and fetch_values commands in MUI that are useful
 *        when just wanting to use MUI to pass data without using spatial
 *        interpolation.
 */

#include "mui.h"
#include "demo8_config.h"

int main( int argc, char ** argv ) {
  if ( argc < 3 ) {
    printf( "USAGE: mpirun -np n1 %s URI1 value1 : -np n2 %s URI2 value2\n\n"
            "n1, n2     : number of ranks for each 'subdomain'\n"
            "URI format : mpi://domain-identifier/interface-identifier\n"
            "value      : an arbitrary number\n\n"
            "EXAMPLE: mpirun -np 1 %s mpi://domain1/ifs 0.618 : -np 1 %s "
            "mpi://domain2/ifs 1.414\n\n",
            argv[0], argv[0], argv[0], argv[0] );
    exit( 0 );
  }

  // Option 1: Declare MUI objects using specialisms (i.e. 1 = 1 dimensional, d = double)
  //mui::uniface1d interface( argv[1] );
  //mui::chrono_sampler_exact1d chrono_sampler;
  //mui::point1d push_point;
  //mui::point1d fetch_point;

  // Option 2: Declare MUI interface and samplers using templates in config.h
  // note: please update types stored in default_config in config.h first to 1-dimensional before compilation
  mui::uniface<mui::demo8_config> interface( argv[1] );
  mui::chrono_sampler_exact<mui::demo8_config> chrono_sampler;
  std::vector<mui::point<mui::demo8_config::REAL, 1>> push_locs(5);

  printf( "domain %s pushed 5 values %s\n", argv[1], argv[2] );

  // Push value stored in "argv[2]" to the MUI interface
  for(size_t i=0; i<5; i++) { //Define push locations as 0-5 and push the value
	  push_locs[i] = static_cast<mui::demo8_config::REAL>(i);
	  interface.push( "data", push_locs[i], atof( argv[2] ) );
  }

  // Commit (transmit by MPI) the values at time=0
  interface.commit( 0 );

  // Fetch the values from the interface using the fetch_points and fetch_values methods
  // (blocking until data at "t=0" exists according to chrono_sampler)
  int time = 0;

  std::vector<mui::point<mui::demo8_config::REAL, 1>> fetch_locs = interface.fetch_points<mui::demo8_config::REAL>( "data", time, chrono_sampler ); // Extract the locations stored in the interface at time=0
  std::vector<double> fetch_vals = interface.fetch_values<mui::demo8_config::REAL>( "data", time, chrono_sampler ); // Extract the values stored in the interface at time=0

  // Print returned values
  for(size_t i=0; i<fetch_locs.size(); i++) {
	  printf( "domain %s fetched value %lf at location %lf\n", argv[1], fetch_vals[i], fetch_locs[i][0] );
  }

  return 0;
}
