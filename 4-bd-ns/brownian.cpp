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

/* flow with a stagnation point in a square box
 * Yazdani, A., Deng, M., Caswell, B., & Karniadakis, G. E. (2016). Flow in complex domains simulated by Dissipative Particle Dynamics driven by geometry-specific body-forces. Journal of Computational Physics, 305, 906-920.
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
 */

int main(int argc, char **argv) {
    using namespace mui;
    MPI_Comm world = mui::mpi_split_by_app();
	uniface2d interface( "mpi://brownian/ifs" );

	// setup parameters
    double box[2] = { -1, 1 };
    double x = 0.1, y = 0.1; // position of the Brownian particle
    double r = 0.1; // particle radius
    double dt = 0.1; // time step size
    double kBT = 1.0; // energy scale
    double eta = 62.9546; // dimensionless viscosity
    double cd  = dt / 6.0 / M_PI / eta / r; // dissipative coefficient
    double cr  = std::sqrt( kBT * dt / 3.0 / M_PI / eta / r ); // fluctuation coefficient
    double drag = 1.0; // coefficient between flow rate and drag force

    // RNG
    std::mt19937 engine( time(0) );
    std::uniform_real_distribution<double> unif( -std::sqrt( 3 * dt ), std::sqrt( 3 * dt ) );

    // trajectory file
    std::ofstream fout( "brownian.txt" );

	// announce first 'span' for smart sending
	interface.announce_recv_span( 0, 1, geometry::sphere2d( {x, y}, r ) );

    // BD run
	for ( int step = 0 ; step <= 10000 ; step++ ) {
        if ( step % 100 == 0 ) printf( "Brownian step %d\n", step );

		// obtain body force exerted by the coupled fluid solver
        sampler_gauss2d<double> s1( r, r/4 );
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

        if ( step % 1 == 0 ) {
        	fout << x << '\t' << y << std::endl;
        }

        // announce updated 'span' for the fluid solver to optimize comm
    	interface.announce_recv_span( step, step+1, geometry::sphere2d( {x, y}, r ) );
        interface.commit( step ); // signaling the other solver that it can move ahead now
    }

    return 0;
}



