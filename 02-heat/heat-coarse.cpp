/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 2                *
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
 * @file heat-coarse.cpp
 * @author Y. H. Tang
 * @date 16 June 2016
 * @brief Coupled simple 1D heat solution with MUI coupling on a coarse
 * grid.
 *
 * Grid scheme, 4:1 coarse-fine ratio, PBC
 * Fine   :                         o-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-o
 *                                  0 1 2 3 4 5 6 7 8 9       ~~~          20
 * Coarse : +-------+-------+-------+-------o-------.-------.-------o-------+-------+-------+-------+
 *          0       1       2       3       4       5       6       7       8       9      10      11
 * +: grid points
 * o: interface points
 * .: place holder for illustrative purposes
 *
 * USAGE: mpirun -np 1 ./heat-coarse : -np 1 ./heat-fine
 */

#include "mui.h"
#include <algorithm>
#include <fstream>

int main( int argc, char ** argv ) {
  const static int N = 12;
  double u1[N], u2[N]; // note that it is not necessary to allocate space for node 5 & 6, but here I do it anyway to
                       // simplify coding
  for ( int i = 0; i <  4; i++ ) u1[i] = 0;
  for ( int i = 8; i < 12; i++ ) u1[i] = 0;

  // Option 1: Declare MUI objects using specialisms (i.e. 1 = 1 dimensional, d = double)
  mui::uniface1d interface( "mpi://coarse/ifs" );
  mui::sampler_gauss1d<double> spatial_sampler( 1, 0.25 );
  mui::temporal_sampler_exact1d temporal_sampler;
  mui::point1d push_point;
  mui::point1d fetch_point;

  // Option 2: Declare MUI objects using templates in config.h
  // note: please update types stored in default_config in config.h first to 1-dimensional before compilation
  //mui::uniface<mui::default_config> interface( "mpi://coarse/ifs" );
  //mui::sampler_gauss<mui::default_config> spatial_sampler( 1, 0.25);
  //mui::temporal_sampler_exact<mui::default_config> temporal_sampler;
  //mui::point<mui::default_config::REAL, 1> push_point;
  //mui::point<mui::default_config::REAL, 1> fetch_point;

  double k = 0.01, H = 1, h = 0.25; // H/h : grid stride for the coarse/fine grid
  double *u = u1, *v = u2;
  std::ofstream fout( "solution-coarse.txt" );

  fout << "TIMESTEP 0" << std::endl;
  for ( int i = 0; i <  4; i++ ) fout << i * H << '\t' << u[i] << '\n';
  for ( int i = 8; i < 12; i++ ) fout << i * H << '\t' << u[i] << '\n';

  for ( int t = 1; t <= 100; t++ ) {
    printf( "Coarse grid step %d\n", t );

    // Push value stored in "u[3]" to the MUI interface
    push_point[0] = 3 * H;
    interface.push( "u", push_point, u[3] );
    // Push value stored in "u[8]" to the MUI interface
    push_point[0] = 8 * H;
    interface.push( "u", push_point, u[8] );
    // Commit (transmit by MPI) the values
    interface.commit( t );

    // Fetch the values for u[4] and u[7] from the interface (blocking until data at "t" exists according to temporal_sampler)
    fetch_point[0] = 4 * H;
    u[4] = interface.fetch( "u", fetch_point, t, spatial_sampler, temporal_sampler );
    fetch_point[0] = 7 * H;
    u[7] = interface.fetch( "u", fetch_point, t, spatial_sampler, temporal_sampler );

    // calculate 'interior' points
    for ( int i = 1; i <  4; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
    for ( int i = 8; i < 11; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
    // calculate 'boundary' points
    v[0]     = u[0    ] + k / ( H * H ) * ( u[1] + u[N - 1] - 2 * u[0    ] );
    v[N - 1] = u[N - 1] + k / ( H * H ) * ( u[0] + u[N - 2] - 2 * u[N - 1] );

    // I/O
    std::swap( u, v );
    fout << "TIMESTEP " << t << std::endl;
    for ( int i = 0; i <  4; i++ ) fout << i * H << '\t' << u[i] << '\n';
    for ( int i = 8; i < 12; i++ ) fout << i * H << '\t' << u[i] << '\n';
  }

  fout.close();

  return 0;
}
