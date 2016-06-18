/*
 * heat-fine.cpp
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

int main( int argc, char ** argv ) {
    using namespace mui;

    const static int N = 21;
	double u1[N], u2[N]; // note that it is not necessary to allocate space for node 0 & 20, but here I do it anyway to simplify coding
    for ( int i = 1; i < 20; i++ ) u1[i] = 2 - i % 2; // spiky I.C.

    uniface1d interface( "mpi://fine/ifs" );

    double k = 0.01, H = 1, h = 0.25; // H/h : grid stride for the coarse/fine grid
    double * u = u1, *v = u2;
    std::ofstream fout( "solution-fine.txt" );

    fout << "TIMESTEP 0" << std::endl;
    for ( int i = 1; i < 20; i++ ) fout << i * h + 3 * H << '\t' << u[i] << '\n';

    for ( int t = 1; t <= 100; t ++ ) {
    	printf("Fine grid step %d\n", t);

    	// push data to the other solver
    	for(int i= 1;i< 8;i++) interface.push( "u", i * h + 3 * H, u[i] );
    	for(int i=13;i<20;i++) interface.push( "u", i * h + 3 * H, u[i] );
    	interface.commit( t );

    	// fetch data from the other solver
    	sampler_exact1d<double> s1;
    	chrono_sampler_exact1d  s2;
    	u[ 0] = interface.fetch( "u",  0 * h + 3 * H, t, s1, s2 );
    	u[20] = interface.fetch( "u", 20 * h + 3 * H, t, s1, s2 );

    	// FDM calculation, all points are 'interior'
    	for(int i=1;i<20;i++) v[i] = u[i] + k / ( h * h ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );

    	// I/O
        std::swap( u, v );
        fout << "TIMESTEP " << t << std::endl;
        for ( int i = 1; i < 20; i++ ) fout << i * h + 3 * H << '\t' << u[i] << '\n';
    }
    fout.close();

    return 0;
}

