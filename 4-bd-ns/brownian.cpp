/*
 * brownian.cpp
 *
 *  Created on: Jul 9, 2016
 *      Author: ytang
 */

#include "../mui/mui.h"
#include <fstream>
#include <limits>
#include <iterator>
#include <cmath>
#include <random>

/* Brownian dynamics in a peridic square box
 * Yazdani, A., Deng, M., Caswell, B., & Karniadakis, G. E. (2016). Flow in complex domains simulated by Dissipative Particle Dynamics driven by geometry-specific body-forces. Journal of Computational Physics, 305, 906-920.
 *                 1
 *     -------------------------
 *    |                         |
 *    |                         |
 *    |                         |
 *    |                         |
 *    |                         |
 *    |                         |
 * -1 |            o            | 1
 *    |                         |
 *    |                         |
 *    |                         |
 *    |                         |
 *    |                         |
 *    |                         |
 *     -------------------------
 *                -1
 *
 * o: intial position of the Brownian particle
 */

int main(int argc, char **argv) {
    using namespace mui;

    MPI_Comm world = mui::mpi_split_by_app();

    double box[2] = { -1, 1 };
    double x = 0.1, y = 0.1; // position of the Brownian particle
    double r = 0.1; // particle radius
    double dt = 0.1; // time step size
    double kBT = 1.0; // energy scale
//    double kBT = 0.0; // energy scale
    double eta = 62.9546; // dimensionless viscosity
    double cd  = dt / 6.0 / M_PI / eta / r; // dissipative coefficient
    double cr  = std::sqrt( kBT * dt / 3.0 / M_PI / eta / r ); // fluctuation coefficient

	uniface2d interface( "mpi://brownian/ifs" );

	std::ofstream fout( "brownian.txt" );

    // RNG
    //std::mt19937 engine( 124 );
    std::mt19937 engine( time(0) );
    std::uniform_real_distribution<double> unif( -std::sqrt( 3 * dt ), std::sqrt( 3 * dt ) );

    // BD run
	geometry::sphere2d recv_region( {x, y}, r );
	interface.announce_recv_span( 0, 1, recv_region );
    for ( int step = 0 ; step <= 10000 ; step++ ) {
        if ( step % 100 == 0 ) printf( "Brownian step %d\n", step );

        double fx = 0, fy = 0;
        double drag = 1.0;

		// body force exerted by the underlying fluid solver
        sampler_gauss2d<double> s1( r, r/4 );
        chrono_sampler_exact2d  s2;
        double ux = interface.fetch( "ux", {x, y}, step, s1, s2 );
        double uy = interface.fetch( "uy", {x, y}, step, s1, s2 );
        //printf("x y %lf %lf ux uy %lf %lf\n", x, y, ux, uy );
        fx = ux * drag;
        fy = uy * drag;

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

//		atoms( i, x ) += ( atoms( i, x ) > box_size[0] ? -box_size[0] : 0.0 ) + ( atoms( i, x ) < 0 ? box_size[0] : 0.0 );
//		if ( atoms( i, y ) < 0    ) atoms( i, y ) = -atoms( i, y );
//		else if ( atoms( i, y ) > box_size[1] ) atoms( i, y ) = 2 * box_size[1] - atoms( i, y );
//		atoms( i, z ) += ( atoms( i, z ) > box_size[2] ? -box_size[2] : 0.0 ) + ( atoms( i, z ) < 0 ? box_size[2] : 0.0 );

        if ( step % 1 == 0 ) {
        	fout << x << '\t' << y << std::endl;
        }

    	geometry::sphere2d recv_region( {x, y}, r );
    	interface.announce_recv_span( step, step+1, recv_region );
        interface.commit( step ); // signaling the other solver that it can move ahead now
    }

    return 0;
}



