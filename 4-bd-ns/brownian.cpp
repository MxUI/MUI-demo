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
 * @file brownian.cpp
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
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <random>

int main( int argc, char ** argv ) {
    using namespace mui;
    MPI_Comm  world = mui::mpi_split_by_app();
    uniface2d interface( "mpi://brownian/ifs" );

    // setup parameters
    double box[2] = {-1, 1};
    double x = 0.1, y = 0.1;                                    // position of the Brownian particle
    double r    = 0.1;                                          // particle radius
    double dt   = 0.1;                                          // time step size
    double kBT  = 1.0;                                          // energy scale
    double eta  = 62.9546;                                      // dimensionless viscosity
    double cd   = dt / 6.0 / M_PI / eta / r;                    // dissipative coefficient
    double cr   = std::sqrt( kBT * dt / 3.0 / M_PI / eta / r ); // fluctuation coefficient
    double drag = 1.0;                                          // coefficient between flow rate and drag force

    // RNG
    std::mt19937                           engine( time( 0 ) );
    std::uniform_real_distribution<double> unif( -std::sqrt( 3 * dt ), std::sqrt( 3 * dt ) );

    // trajectory file
    std::ofstream fout( "brownian.txt" );

    // announce first 'span' for smart sending
    interface.announce_recv_span( 0, 1, geometry::sphere2d( {x, y}, r ) );

    // BD run
    for ( int step = 0; step <= 10000; step++ ) {
        if ( step % 100 == 0 ) printf( "Brownian step %d\n", step );

        // obtain body force exerted by the coupled fluid solver
        sampler_gauss2d<double> s1( r, r / 4 );
        chrono_sampler_exact2d  s2;
        auto ux = interface.fetch( "ux", {x, y}, step, s1, s2 );
        auto uy = interface.fetch( "uy", {x, y}, step, s1, s2 );
        auto fx = ux * drag;
        auto fy = uy * drag;

        // BD overdamped integration & temperature calculation
        auto Bx = unif( engine );
        auto By = unif( engine );
        auto dx = cd * fx + cr * Bx;
        auto dy = cd * fy + cr * By;
        x += dx;
        y += dy;

        // an empty line tells gnuplot to not connect the adjacent points
		if ( x > box[1] ) { x -= box[1] - box[0]; fout << std::endl; }
		if ( x < box[0] ) { x += box[1] - box[0]; fout << std::endl; }
		if ( y > box[1] ) { y -= box[1] - box[0]; fout << std::endl; }
		if ( y < box[0] ) { y += box[1] - box[0]; fout << std::endl; }

        if ( step % 1 == 0 ) { fout << x << '\t' << y << std::endl; }

        // announce updated 'span' for the fluid solver to optimize comm
        interface.announce_recv_span( step, step + 1, geometry::sphere2d( {x, y}, r ) );
        interface.commit( step ); // signaling the other solver that it can move ahead now
    }

    return 0;
}
