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
 * @file sph-standalone.cpp
 * @author X. Bian
 * @date 19 June 2016
 * @brief Stand alone SPH solver for hybrid Lagrangian-Eulerian simulation
 * of Heat Conduction.
 */

#include <cstdio>

double cubic_spline_gradient( double r, double rc ) {
    const static double sigma  = 2.0 / 3.0; // for one diemension
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

    const int    N = 19, Nint = 3;      // number of inside/outside particles
    const int    Ntotal = N + Nint * 2; // number of total particles
    const double L      = 10;           // box length for internal materia

    double u[Ntotal], x[Ntotal], du[Ntotal];
    double rho   = 1000.0; // density and mass: constants
    double kappa = 1.0;    // conductivity: constant
    double cv    = 1.0;    // heat compacity at consant volume
    double dx = L / ( N + 1 ), rc = 2.5 * dx;
    double m      = dx * rho;
    double T_left = 0.0, T_right = 1.0;
    double t_end = 0.5, dt = 0.25 * dx * dx / kappa;
    int    step_end = 1000;

    // initial conditions
    for ( int i = 0; i < Ntotal; i++ ) {
        x[i] = 0.0 + ( i + 1 - Nint ) * dx;
        u[i] = cv * ( i <= Ntotal / 2 ? T_left : T_right );
    }

    for ( int i = 0; i < Ntotal; i++ ) printf( "%f, %f\n", x[i], u[i] );
    printf( "\n" );

    // time integration forward Euler scheme
    for ( int k = 0; k < step_end; k++ ) {

        // reset du
        for ( int i = 0; i < Ntotal; i++ ) { du[i] = 0.0; }

        // N^2 brute-force pairwise evaluation
        for ( int i = 0; i < Ntotal - 1; i++ ) {
            for ( int j = i + 1; j < Ntotal; j++ ) {
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
        for ( int i = Nint; i < Ntotal - Nint; i++ ) u[i] += du[i];

        if ( k % 100 == 0 ) {
            for ( int i = 0; i < Ntotal; i++ ) printf( "%f, %f\n", x[i], u[i] );
            printf( "\n" );
        }
    }

    return 0;
}
