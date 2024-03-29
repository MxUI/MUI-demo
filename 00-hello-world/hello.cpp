/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 0                *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
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
 * @file hello.cpp
 * @author Y. H. Tang
 * @date 9 June 2016
 * @brief Simple "Hello World" example to demonstrate MUI.
 */

#include "mui.h"

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
  mui::uniface1d interface( argv[1] );
  mui::sampler_exact1d<double> spatial_sampler;
  mui::temporal_sampler_exact1d temporal_sampler;
  mui::point1d push_point;
  mui::point1d fetch_point;

  // Option 2: Declare MUI interface and samplers using templates in config.h
  // note: please update types stored in default_config in config.h first to 1-dimensional before compilation
  //mui::uniface<mui::default_config> interface( argv[1] );
  //mui::sampler_exact<mui::default_config> spatial_sampler;
  //mui::temporal_sampler_exact<mui::default_config> temporal_sampler;
  //mui::point<mui::default_config::REAL, 1> push_point;
  //mui::point<mui::default_config::REAL, 1> fetch_point;

  printf( "domain %s pushed value %s\n", argv[1], argv[2] );

  // Push value stored in "argv[2]" to the MUI interface
  push_point[0] = 0;
  interface.push( "data", push_point, atof( argv[2] ) );

  // Commit (transmit by MPI) the value
  interface.commit( 0 );

  // Fetch the value from the interface (blocking until data at "t=0" exists according to temporal_sampler)
  int time = 0;
  fetch_point[0] = 0;
  double v = interface.fetch( "data", fetch_point, time, spatial_sampler, temporal_sampler );

  printf( "domain %s fetched value %lf\n", argv[1], v );

  return 0;
}
