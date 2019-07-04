/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 3                *
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
 * @file fdm.cpp
 * @author Y. H. Tang
 * @date 16 June 2016
 * @brief Hybrid Lagrangian-Eulerian Simulation of Heat Conduction.
 *
 * SPH-Finite Difference, PBC
 * SPH :                         o o o o o o o o o o o o o o o o o o o o o
 *                               0 1 2 3 4 5 6 7 8 9      .....         24
 *                               ^ ^   * *                       * *   ^ ^
 * FDM :    +---+---+---+---+---+---+---+                         +---+---+---+---+---+---+---+
 *          0   1   2   3   4   5   6   7                         8   9  10  11  12  13  14  15
 *                              *   *   ^                         ^   *   *
 * +: grid points
 * o: SPH particles
 * ^: interface points - fetch
 * *: interface points - push
 *
 * // USAGE: mpirun -np 1 ./sph : -np 1 ./fdm
 */

#include "../mui/mui.h"
#include <algorithm>
#include <fstream>

int main( int argc, char ** argv ) {

    using namespace mui;
    uniface1d interface( "mpi://fdm/ifs" );

    const static int N = 16;
    double           k = 0.0625, H = 1; // H: grid stride for the coarse/fine grid
    auto i2x = [&]( int i ) {           // convert grid index to position
        if ( i < 8 )
            return -11.25 + i * H;
        else
            return 4.25 + ( i - 8 ) * H;
    };
    double u1[N], u2[N]; // note that it is not necessary to allocate space for node 5 & 6, but here I do it anyway to
                         // simplify coding
    for ( int i = 0; i < 8; i++ ) u1[i] = 0;
    for ( int i = 8; i < 16; i++ ) u1[i] = 1;
    double *  u = u1, *v = u2;

    std::ofstream fout( "solution-fdm.txt" );

    for ( int t = 0; t < 100; t++ ) {
        // push data to the other solver
        for ( int i : {5, 6, 9, 10} ) interface.push( "u", i2x( i ), u[i] );
        interface.commit( t );

        // fetch data from the other solver
        sampler_moving_average1d<double> avg( 1 );
        chrono_sampler_exact1d           exact;
        for ( int i : {7, 8} ) u[i] = interface.fetch( "u", i2x( i ), t, avg, exact );

        // calculate 'interior' points
        for ( int i = 1; i < 7; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        for ( int i = 9; i < 16; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        // calculate 'boundary' points
        v[0]     = u[0] + k / ( H * H ) * ( u[1] + u[N - 1] - 2 * u[0] );
        v[N - 1] = u[N - 1] + k / ( H * H ) * ( u[0] + u[N - 2] - 2 * u[N - 1] );

        // I/O
        std::swap( u, v );
        if ( t % 10 == 0 ) {
            printf( "FDM step %d\n", t );
            for ( int i = 0; i < 7; i++ ) fout << i2x( i ) << '\t' << u[i] << '\n';
            for ( int i = 9; i < 16; i++ ) fout << i2x( i ) << '\t' << u[i] << '\n';
            fout << std::endl;
        }
    }
    fout.close();

    return 0;
}
