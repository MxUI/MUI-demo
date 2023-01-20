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
 * @file heat-fine.cpp
 * @author Y. H. Tang
 * @date 16 June 2016
 * @brief Coupled simple 1D heat solution with MUI coupling on a fine
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
  const static int N = 21;
  double u1[N], u2[N];

  for ( int i = 1; i < 20; i++ ) u1[i] = 2 - i % 2 + ( i - 1 ) * 0.1; // spiky, skewed I.C.

  // Option 1: Declare MUI objects using specialisms (i.e. 1 = 1 dimensional, d = double)
  mui::uniface1d interface( "mpi://fine/ifs" );
  mui::sampler_exact1d<double> spatial_sampler;
  mui::chrono_sampler_exact1d chrono_sampler;
  mui::point1d push_point;
  mui::point1d fetch_point;

  // Option 2: Declare MUI objects using templates in config.h
  // note: please update types stored in default_config in config.h first to 1-dimensional before compilation
  //mui::uniface<mui::default_config> interface( "mpi://fine/ifs" );
  //mui::sampler_exact<mui::default_config> spatial_sampler;
  //mui::chrono_sampler_exact<mui::default_config> chrono_sampler;
  //mui::point<mui::default_config::REAL, 1> push_point;
  //mui::point<mui::default_config::REAL, 1> fetch_point;

  double        k = 0.01, H = 1, h = 0.25; // H/h : grid stride for the coarse/fine grid
  double *      u = u1, *v = u2;
  std::ofstream fout( "solution-fine.txt" );

  fout << "TIMESTEP 0" << std::endl;
  for ( int i = 1; i < 20; i++ ) fout << i * h + 3 * H << '\t' << u[i] << '\n';

  for ( int t = 1; t <= 100; t++ ) {
    printf( "Fine grid step %d\n", t );

    // Push values to the MUI interface
    for ( int i =  1; i <  8; i++ ) {
      push_point[0] = i * h + 3 * H;
      interface.push( "u", push_point, u[i] );
    }
    for ( int i = 13; i < 20; i++ ) {
      push_point[0] = i * h + 3 * H;
      interface.push( "u", push_point, u[i] );
    }
    // Commit (transmit by MPI) the values
    interface.commit( t );

    // Fetch the values from the interface (blocking until data at "t" exists according to chrono_sampler)
    fetch_point[0] = 0 * h + 3 * H;
    u[0] = interface.fetch( "u",  fetch_point, t, spatial_sampler, chrono_sampler );
    fetch_point[0] = 20 * h + 3 * H;
    u[20] = interface.fetch( "u", fetch_point, t, spatial_sampler, chrono_sampler );

    // FDM calculation, all points are 'interior'
    for ( int i = 1; i < 20; i++ ) v[i] = u[i] + k / ( h * h ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );

    // I/O
    std::swap( u, v );
    fout << "TIMESTEP " << t << std::endl;
    for ( int i = 1; i < 20; i++ ) fout << i * h + 3 * H << '\t' << u[i] << '\n';
  }

  fout.close();

  return 0;
}
