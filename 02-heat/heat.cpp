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
 * @file heat.cpp
 * @author Y. H. Tang
 * @date 16 June 2016
 * @brief Standalone reference simple 1D heat solution without MUI coupling.
 */

#include <algorithm>
#include <fstream>

int main( int argc, char ** argv ) {
    const static int N = 16;
    double           u1[N], u2[N];
    for ( int i = 0; i < N; i++ ) u1[i] = ( i >= N / 2 ) ? 1 : 0;

    double        k = 0.1, h = 1;
    double *      u = u1, *v = u2;
    std::ofstream fout( "solution.txt" );

    fout << "TIMESTEP 0" << std::endl;
    for ( int i = 0; i < N; i++ ) fout << i * h << '\t' << u[i] << '\n';

    for ( int t = 1; t <= 100; t++ ) {
        for ( int i = 1; i < N - 1; i++ ) v[i] = u[i] + k / ( h * h ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        v[0]                                   = u[0] + k / ( h * h ) * ( u[1] + u[N - 1] - 2 * u[0] );
        v[N - 1]                               = u[N - 1] + k / ( h * h ) * ( u[0] + u[N - 2] - 2 * u[N - 1] );

        std::swap( u, v );
        fout << "TIMESTEP " << t << std::endl;
        for ( int i = 0; i < N; i++ ) fout << i * h << '\t' << u[i] << '\n';
    }
    fout.close();

    return 0;
}
