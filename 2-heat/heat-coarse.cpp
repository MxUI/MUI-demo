/*
 * heat-coarse.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: ytang
 */


#include "../mui/mui.h"
#include <algorithm>
#include <fstream>

/* Grid scheme, 4:1 coarse-fine ratio, PBC
 * Fine   :                         o-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-o
 *                                  0 1 2 3 4 5 6 7 8 9       ~~~          20
 * Coarse : +-------+-------+-------+-------o-------.-------.-------o-------+-------+-------+-------+
 *          0       1       2       3       4       5       6       7       8       9      10      11
 * +: grid points
 * o: interface points
 * .: place holder for illustrative purposes
 */

// USAGE: mpirun -np 1 ./heat-coarse : -np 1 ./heat-fine

int main( int argc, char ** argv ) {
    using namespace mui;

    const static int N = 12;
    double u1[N], u2[N]; // note that it is not necessary to allocate space for node 5 & 6, but here I do it anyway to simplify coding
    for ( int i = 0; i <  4; i++ ) u1[i] = 0;
    for ( int i = 8; i < 12; i++ ) u1[i] = 0;

    uniface1d interface( "mpi://coarse/ifs" );

    double k = 0.01, H = 1, h = 0.25; // H/h : grid stride for the coarse/fine grid
    double * u = u1, *v = u2;
    std::ofstream fout( "solution-coarse.txt" );

    fout << "TIMESTEP 0" << std::endl;
    for ( int i = 0; i <  4; i++ ) fout << i * H << '\t' << u[i] << '\n';
    for ( int i = 8; i < 12; i++ ) fout << i * H << '\t' << u[i] << '\n';

    for ( int t = 1; t <= 100; t ++ ) {
        printf( "Coarse grid step %d\n", t );

        // push data to the other solver
        interface.push( "u", 3 * H, u[3] );
        interface.push( "u", 8 * H, u[8] );
        interface.commit( t );

        // fetch data from the other solver
        sampler_gauss1d<double> s1( 1, 0.25 );
        chrono_sampler_exact1d  s2;
        u[4] = interface.fetch( "u", 4 * H, t, s1, s2 );
        u[7] = interface.fetch( "u", 7 * H, t, s1, s2 );

        // calculate 'interior' points
        for ( int i = 1; i < 4; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        for ( int i = 8; i < 11; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        // calculate 'boundary' points
        v[    0] = u[    0] + k / ( H * H ) * ( u[1] + u[N - 1] - 2 * u[    0] );
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

