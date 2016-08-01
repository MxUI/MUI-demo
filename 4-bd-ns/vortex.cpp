/*
 * vortex.cpp
 *
 *  Created on: Jul 9, 2016
 *      Author: ytang
 */

#include "../mui/mui.h"
#include <algorithm>
#include <fstream>

// USAGE: mpirun -np 1 ./brownian : -np 4 ./vortex

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

int main( int argc, char ** argv ) {
    using namespace mui;
    MPI_Comm world = mui::mpi_split_by_app();
	uniface2d interface( "mpi://vortex/ifs" );

    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );
    if ( size != 4 ) { // hard-coded the box to be decomposed into 2-by-2 ranks
    	printf("This solver must be ran with 4 MPI ranks\n");
    	exit(0); // no need to MPI_Finalize - MUI hook will take care
    }
    int rank_x = rank / 2;
    int rank_y = rank % 2;

    // decompose the domain
    constexpr static int N = 40; // number of grid points in each local domain
	constexpr static double h = 1.0 / N;
	double local_x0 = rank_x * 1 + -1; // local origin
	double local_y0 = rank_y * 1 + -1;
	double local_x1 = local_x0 + 1;
	double local_y1 = local_y0 + 1;
	double px[N][N], py[N][N], ux[N][N], uy[N][N];

	// generate 'fake' flow field from analytical solution
    double Gamma = 100; // vortex circulation
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++) {
			double x = i * h + h/2 + local_x0;
			double y = j * h + h/2 + local_y0;
			px[i][j] = x;
			py[i][j] = y;
			ux[i][j] = 0;
			uy[i][j] = 0;
			for(auto bx : {-1,1})
				for( auto by: {-1,1} ) {
					auto dx = x - bx;
					auto dy = y - by;
					auto r = std::sqrt( dx * dx + dy * dy );
					auto u = Gamma / 2.0 / M_PI / r;
					ux[i][j] += u * -dy / r;
					uy[i][j] += u *  dx / r;
				}
		}

	// dump the flow field
	std::stringstream out;
	for(int r=0;r<4;r++) {
		if ( rank == r ) {
			printf("rank %d\n", rank);
			for(int i=0;i<N;i++)
				for(int j=0;j<N;j++)
					out << px[i][j] << '\t' << py[i][j] << '\t' << ux[i][j] << '\t' << uy[i][j] << std::endl;
		}
	}
	MPI_File fout;
	MPI_File_open( world, "vortex.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fout );
	MPI_Offset bytes, offset;
	bytes = out.str().size();
	MPI_Exscan( &bytes, &offset, 1, MPI_OFFSET, MPI_SUM, world );
	MPI_Status status;
	printf("rank %d size %d offset %d\n", rank, bytes, offset);
	MPI_File_write_at_all( fout, offset, out.str().c_str(), bytes, MPI_CHAR, &status );
	MPI_File_close( &fout );

	// annouce send span
	geometry::box2d send_region( {local_x0, local_y0}, {local_x1, local_y1} );
	printf("send region for rank %d: %lf %lf - %lf %lf\n", rank, local_x0, local_y0, local_x1, local_y1);
	interface.announce_send_span( 0, 10000, send_region );

    for ( int t = 0; t <= 10000; t ++ ) {
        // push data to the other solver
        for(int i=0;i<N;i++) {
        		for(int j=0;j<N;j++) {
        			point2d loc( px[i][j], py[i][j] );
        			interface.push( "ux", loc, ux[i][j] );
        			interface.push( "uy", loc, uy[i][j] );
        		}
        }

        int sent = interface.commit( t );
        if ( t % 100 == 0 ) printf("Vortex rank %d step %d actual sending: %s\n", rank, t, sent ? "ON" : "OFF" );

        // no need to fetch data, but wait for the other side to catch up
        interface.barrier( t );
    }

    return 0;
}




