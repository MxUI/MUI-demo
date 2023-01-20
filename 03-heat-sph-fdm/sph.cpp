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
 * @file sph.cpp
 * @author X. Bian & Y. H. Tang
 * @date 19 June 2016
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
 * USAGE: mpirun -np 1 ./sph : -np 1 ./fdm
 */

#include "mui.h"
#include <cstdio>
#include <fstream>

double cubic_spline_gradient( double r, double rc ) {
    const static double sigma  = 2.0 / 3.0; // for one dimensional problem
    double              h      = rc / 2.0;
    double              coef_g = sigma / h / h;
    double              s      = r / h;
    double              w_g    = 0;
    if ( s < 1.0 )
        w_g = -3.0 * s + 9.0 / 4.0 * s * s;
    else if ( s < 2.0 )
        w_g = -3.0 / 4.0 * ( 2 - s ) * ( 2 - s );
    return w_g * coef_g;
}

int main() {
    const int    Ni = 21, No = 2; // number of inside/outside particles
    const int    N = Ni + No * 2; // number of total particles
    const double L = 11;          // box length for internal materia

    double u[N], x[N], du[N];
    double rho   = 100.0; // density and mass: constants
    double kappa = 1.0;   // conductivity: constant
    double cv    = 1.0;   // heat compacity at consant volume
    double dx = L / ( Ni + 1 ), rc = 2.5 * dx;
    double m  = dx * rho;
    double dt = 0.25 * dx * dx / kappa;

    // Option 1: Declare MUI objects using specialisms (i.e. 1 = 1 dimensional, d = double)
    mui::uniface1d interface( "mpi://sph/ifs" );
    mui::sampler_gauss1d<double> spatial_sampler( rc, rc / 2 );
    mui::chrono_sampler_exact1d chrono_sampler;
    mui::point1d push_point;
    mui::point1d fetch_point;

    // Option 2: Declare MUI objects using templates in config.h
    // note: please update types stored in default_config in config.h first to 1-dimensional before compilation
    //mui::uniface<mui::default_config> interface( "mpi://sph/ifs" );
    //mui::sampler_gauss<mui::default_config> spatial_sampler( rc, rc / 2);
    //mui::chrono_sampler_exact<mui::default_config> chrono_sampler;
    //mui::point<mui::default_config::REAL, 3> push_point;
    //mui::point<mui::default_config::REAL, 3> fetch_point;

    // initial conditions
    for ( int i = 0; i < N; i++ ) {
      x[i] = ( i - N / 2 ) * dx;
      u[i] = cv * ( i <= N / 2 ? 0 : 1 );
      if ( i >= N / 4 && i <= N * 3 / 4 ) u[i] += i % 2 * 0.5;
    }

    std::ofstream fout( "solution-sph.txt" );

    // time integration forward Euler scheme
    for ( int k = 0; k < 100; k++ ) {
      // Push values to the MUI interface
      for ( int i : {3, 4, N - 5, N - 4} ) {
        push_point[0] = x[i];
        interface.push( "u", push_point, u[i] );
      }
      // Commit (transmit by MPI) the values
      interface.commit( k );

      // Fetch the values from the interface (blocking until data at "k" exists according to chrono_sampler)
      for ( int i : {0, 1, N - 2, N - 1} ) {
        fetch_point[0] = x[i];
        u[i] = interface.fetch( "u", fetch_point, k, spatial_sampler, chrono_sampler );
      }

      // reset du
      for ( int i = 0; i < N; i++ ) du[i] = 0.0;

      // N^2 brute-force pairwise evaluation
      for ( int i = 0; i < N - 1; i++ ) {
        for ( int j = i + 1; j < N; j++ ) {
          double r_ij = x[j] - x[i];

          if ( r_ij <= rc ) {
            double w_g   = cubic_spline_gradient( r_ij, rc );
            double du_ij = 2.0 * kappa * m / rho / rho * w_g * ( u[i] / cv - u[j] / cv );
            du[i] += du_ij;
            du[j] -= du_ij;
          }
        }
      }

      // only update the internal material points
      for ( int i = No; i < N - No; i++ ) u[i] += du[i];

      if ( k % 10 == 0 ) {
        printf( "SPH step %d\n", k );
        for ( int i = No; i < N - No; i++ ) fout << x[i] << '\t' << u[i] << std::endl;
        fout << std::endl;
      }
    }

    return 0;
}
