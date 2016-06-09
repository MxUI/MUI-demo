/*
 * hello.cpp
 *
 *  Created on: Jun 9, 2016
 *      Author: ytang
 */

#include "../mui/mui.h"

int main( int argc, char ** argv ) {
    using namespace mui;

    if ( argc < 3 ) {
        printf( "USAGE: mpirun -np n1 %s URI1 value1 : -np n2 %s URI2 value2\n\n"
                "n1, n2     : number of ranks for each 'subdomain'\n"
                "URI format : mpi://domain-identifier/interface-identifier\n"
                "value      : an arbitrary number\n\n"
                "EXAMPLE: mpirun -np 1 %s mpi://domain1/ifs 0.618 : -np 1 %s mpi://domain2/ifs 1.414\n\n",
                argv[0], argv[0], argv[0], argv[0] );
        exit( 0 );
    }

    uniface1d interface( argv[1] );
    printf( "domain %s pushed value %s\n", argv[1], argv[2] );
    interface.push( "data", 0, atof( argv[2] ) );
    interface.commit( 0 );
    double v = interface.fetch( "data", 0, 0, sampler_exact1d<double>(), chrono_sampler_exact1d() );
    printf( "domain %s fetched value %lf\n", argv[1], v );

    return 0;
}


