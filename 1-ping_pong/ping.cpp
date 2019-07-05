/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 1                *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
* R. Sawko*                                                                  *
* (*IBM Research)                                                            *
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
 * @file ping.cpp
 * @author Robert Sawko
 * @date 10 June 2016
 * @brief Classical ping-pong MPI communication demonstration using MUI as
 * the interface.
 *
 * USAGE: mpirun -np 1 ./ping : -np 1 ./pong
 */

#include "../mui/mui.h"

int main() {
  // Option 1: Declare MUI objects using specialisms (i.e. 1 = 1 dimensional, d = double)
  mui::uniface1d interface( "mpi://ping/ifs" );
  mui::sampler_exact1d<int> spatial_sampler;
  mui::chrono_sampler_exact1d chrono_sampler;
  mui::point1d push_point;
  mui::point1d fetch_point;

  // Option 2: Declare MUI objects using templates in config.h
  // note: please update types stored in default_config in config.h first to 1-dimensional before compilation
  //mui::uniface<mui::default_config> interface( "mpi://ping/ifs" );
  //mui::sampler_exact<mui::default_config> spatial_sampler;
  //mui::chrono_sampler_exact<mui::default_config> chrono_sampler;
  //mui::point<mui::default_config::REAL, 1> push_point;
  //mui::point<mui::default_config::REAL, 1> fetch_point;

  int state = 0;

  for ( int t = 0; t < 100; ++t ) {
    state++;
    // Push value stored in "state" to the MUI interface
    push_point[0] = 0;
    interface.push( "data", push_point, state );
    // Commit (transmit by MPI) the value
    interface.commit( t );
    // Fetch the value from the interface (blocking until data at "t" exists according to chrono_sampler)
    fetch_point[0] = 0;
    state = interface.fetch( "data", fetch_point, t, spatial_sampler, chrono_sampler );
  }

  printf( "Final ping state: %d\n", state );

  return 0;
}
