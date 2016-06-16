/*
 * heat.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: ytang
 */


//#include "../mui/mui.h"
#include <algorithm>
#include <fstream>

int main( int argc, char ** argv ) {
    const static int N = 16;
    double u1[N], u2[N];
    for ( int i = 0; i < N; i++ ) u1[i] = ( i >= N / 2 ) ? 1 : 0;

    double k = 0.1, h = 1;
    double * u = u1, *v = u2;
    std::ofstream fout( "solution.txt" );
    fout << "TIMESTEP 0" << std::endl;
    for ( int i = 0; i < N; i++ ) fout << i * h << '\t' << u[i] << '\n';
    for ( int t = 1; t <= 100; t ++ ) {
        for ( int i = 1; i < N - 1; i++ ) v[i] = u[i] + k / ( h * h ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        v[  0] = u[  0] + k / ( h * h ) * ( u[1] + u[N - 1] - 2 * u[  0] );
        v[N - 1] = u[N - 1] + k / ( h * h ) * ( u[0] + u[N - 2] - 2 * u[N - 1] );
        std::swap( u, v );
        fout << "TIMESTEP " << t << std::endl;
        for ( int i = 0; i < N; i++ ) fout << i * h << '\t' << u[i] << '\n';
    }
    fout.close();

    return 0;
}

