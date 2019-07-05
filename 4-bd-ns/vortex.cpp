/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 4                *
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
 * @file vortex.cpp
 * @author Y. H. Tang
 * @date 9 July 2016
 * @brief Flow with a stagnation point in a square box.
 *
 * Yazdani, A., Deng, M., Caswell, B., & Karniadakis, G. E. (2016).
 * Flow in complex domains simulated by Dissipative Particle Dynamics
 * driven by geometry-specific body-forces. Journal of Computational Physics,
 * 305, 906-920.
 *
 *                 1
 *     -------------------------
 *    |      |    | |    |      |
 *    |      .    . .    .      |
 *    |     .    .   .    .     |
 *    |-<--'    .     .    '-->-|
 *    |        .       .        |
 *    |---<---'         '--->---|
 * -1 |            X o          | 1
 *    |---<---.         .--->---|
 *    |        '       '        |
 *    |-<--.    '     '    .-->-|
 *    |     '    '   '    '     |
 *    |      '    ' '    '      |
 *    |       |   | |    |      |
 *     -------------------------
 *                -1
 *
 * X: stagnation point
 * o: particle initial position
 *
 * USAGE: mpirun -np 1 ./brownian : -np 4 ./vortex
 */

#include "../mui/mui.h"
#include <algorithm>
#include <fstream>

int main( int argc, char ** argv ) {
  MPI_Comm  world = mui::mpi_split_by_app();

  // Option 1: Declare MUI objects using specialisms (i.e. 2 = 2 dimensional, d = double)
  mui::uniface2d interface( "mpi://vortex/ifs" );
  mui::point2d push_point;

  // Option 2: Declare MUI interface and samplers using templates in config.h
  // note: please update types stored in default_config in config.h first to 2-dimensional before compilation
  //mui::uniface<mui::default_config> interface( "mpi://brownian/ifs" );
  //mui::point<mui::default_config::REAL, 2> push_point;

  int rank, size;
  MPI_Comm_rank( world, &rank );
  MPI_Comm_size( world, &size );
  if ( size != 4 ) { // hard-coded the box to be decomposed into 2-by-2 ranks
    printf( "This solver must be ran with 4 MPI ranks\n" );
    exit( 0 ); // no need to MPI_Finalize - MUI hook will take care
  }
  int rank_x = rank / 2;
  int rank_y = rank % 2;

  // decompose the domain
  constexpr static int N = 40; // number of grid points in each local domain
  constexpr static double h = 1.0 / N;
  double local_x0 = rank_x * 1 + -1; // local origin
  double local_y0 = rank_y * 1 + -1;
  double local_x1 = local_x0 + 1;
  double local_y1 = local_y0 + 1;
  double px[N][N], py[N][N], ux[N][N], uy[N][N];

  // generate 'fake' flow field from analytical solution
  double Gamma = 100; // vortex circulation
  for ( int i = 0; i < N; i++ ) {
    for ( int j = 0; j < N; j++ ) {
      double x = i * h + h / 2 + local_x0;
      double y = j * h + h / 2 + local_y0;
      px[i][j] = x;
      py[i][j] = y;
      ux[i][j] = 0;
      uy[i][j] = 0;
      for ( auto bx : {-1, 1} ) {
        for ( auto by : {-1, 1} ) {
          auto dx = x - bx;
          auto dy = y - by;
          auto r  = std::sqrt( dx * dx + dy * dy );
          auto u  = Gamma / 2.0 / M_PI / r;
          ux[i][j] += u * -dy / r;
          uy[i][j] += u * dx / r;
        }
      }
    }
  }

  // dump the flow field
  std::stringstream out;
  for ( int r = 0; r < 4; r++ ) {
    if ( rank == r ) {
      printf( "rank %d\n", rank );
      for ( int i = 0; i < N; i++ ) {
        for ( int j = 0; j < N; j++ ) {
          out << px[i][j] << '\t' << py[i][j] << '\t' << ux[i][j] << '\t' << uy[i][j] << std::endl;
        }
      }
    }
  }

  MPI_File fout;
  MPI_File_open( world, "vortex.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fout );
  MPI_Offset bytes, offset;
  bytes = out.str().size();
  MPI_Exscan( &bytes, &offset, 1, MPI_OFFSET, MPI_SUM, world );
  MPI_Status status;
  printf( "rank %d size %d offset %d\n", rank, bytes, offset );
  MPI_File_write_at_all( fout, offset, out.str().c_str(), bytes, MPI_CHAR, &status );
  MPI_File_close( &fout );

  // annouce send span
  mui::geometry::box2d send_region( {local_x0, local_y0}, {local_x1, local_y1} );
  printf( "send region for rank %d: %lf %lf - %lf %lf\n", rank, local_x0, local_y0, local_x1, local_y1 );
  interface.announce_send_span( 0, 10000, send_region );

  for ( int t = 0; t <= 10000; t++ ) {
    // push data to the other solver
    for ( int i = 0; i < N; i++ ) {
      for ( int j = 0; j < N; j++ ) {
        push_point[0] = px[i][j];
        push_point[1] = py[i][j];
        interface.push( "ux", push_point, ux[i][j] );
        interface.push( "uy", push_point, uy[i][j] );
      }
    }

    int sent = interface.commit( t );

    if ( t % 100 == 0 ) printf( "Vortex rank %d step %d actual sending: %s\n", rank, t, sent ? "ON" : "OFF" );

    // no need to fetch data, but wait for the other side to catch up
    interface.barrier( t );
  }

  return 0;
}
