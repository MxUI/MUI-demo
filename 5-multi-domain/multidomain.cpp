/*
 * multidomain.cpp
 *
 *  Created on: Oct 18, 2016
 *      Author: ytang
 */

#include "../mui/mui.h"

/* Test establishing multiple uniface instances using several inter-communicators
 *
 * USAGE 1: make
 * USAGE 2: after compilation, run:
 *   mpirun -np n1 ./multidomain domain_name1 interface1 interface2 ... : \
 *          -np n2 ./multidomain domain_name2 interface2 interface3 ... : \
 *          ...
 */

int main( int argc, char ** argv ) {
    mui::mpi_split_by_app();

    std::string              domain = argv[1];
    std::vector<std::string> interfaces;
    for ( int i = 2; i < argc; i++ ) interfaces.emplace_back( argv[i] );

    auto ifs = mui::create_uniface<mui::config_3d>( domain, interfaces );
}
