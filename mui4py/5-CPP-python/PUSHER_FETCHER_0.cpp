/**
 *
 * @file Pusher_Fetcher_0.cpp
 * @author W.L.
 * @date 10 Jan 2019
 * @brief Pusher_Fetcher_0.cpp <-> MUI (C++) <-> MUI (Python) <-> Pusher_Fetcher_1.py two way Coupling test case.
 *
 */

#include "mui.h"

int main(int argc, char ** argv) {
    using namespace mui;

	uniface3d interface( "mpi://PUSHER_FETCHER_0/COUPLING"  );

	MPI_Comm  world = mui::mpi_split_by_app();

    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

	// setup parameters
    constexpr static int    Nx        = 1; // number of grid points in x axis
    constexpr static int    Ny        = 11; // number of grid points in y axis
    constexpr static int    Nz        = 1; // number of grid points in z axis
	const char* name_push = "data_cpp";		
	const char* name_fetch = "data_python";
    double r    = 1e-10;                      // search radius	
    int Nt = Nx * Ny * Nz; // total time steps	
    int steps = 10; // total time steps
	double local_x0 = 6; // local origin
    double local_y0 = 0;
	double local_z0 = 0;
    double local_x1 = 6;
    double local_y1 = 1;
	double local_z1 = 0;
	double local_x2 = 4; // local origin
    double local_y2 = 0;
	double local_z2 = 0;
    double local_x3 = 4;
    double local_y3 = 1;
	double local_z3 = 0;
    double pp[Nx][Ny][Nz][3], pf[Nx][Ny][Nz][3], cpp_push[Nx][Ny][Nz], python_fetch[Nx][Ny][Nz];

	// Push points generation and evaluation
	for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int k = 0; k < Nz; ++k ) {
				double x = 6.0;
				double y = j * 0.1;
				double z = 0.0;
				pp[i][j][k][0] = x;
				pp[i][j][k][1] = y;
				pp[i][j][k][2] = z;
				cpp_push[i][j][k] = 32.0;
			}
        }
	}

	// Fetch points generation and evaluation
	for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int k = 0; k < Nz; ++k ) {
				double x = 4.0;
				double y = j * 0.1;
				double z = 0.0;
				pf[i][j][k][0] = x;
				pf[i][j][k][1] = y;
				pf[i][j][k][2] = z;
				python_fetch[i][j][k] = 0.0;
			}
        }
	}

    // annouce send span
    geometry::box3d send_region( {local_x0, local_y0, local_z0}, {local_x1, local_y1, local_z1} );
    geometry::box3d recv_region( {local_x2, local_y2, local_z2}, {local_x3, local_y3, local_z3} );
    printf( "{PUSHER_FETCHER_0} send region for rank %d: %lf %lf %lf - %lf %lf %lf\n", rank, local_x0, local_y0, local_z0, local_x1, local_y1, local_z1 );
    interface.announce_send_span( 0, steps, send_region );
    interface.announce_recv_span( 0, steps, recv_region );

	// define spatial and temporal samplers
	sampler_gauss3d<double> s1( r, r / 4 );
	chrono_sampler_exact3d  s2;
	
	// commit ZERO step
	interface.commit(0);

	// Begin time loops
    for ( int n = 0; n < steps; ++n ) {

		printf("\n");
		printf("{PUSHER_FETCHER_0} %d Step ", n );	

		// push data to the other solver
		for ( int i = 0; i < Nx; ++i ) {
            for ( int j = 0; j < Ny; ++j ) {
				for ( int k = 0; k < Nz; ++k ) {
					point3d locp( pp[i][j][k][0], pp[i][j][k][1], pp[i][j][k][2] );
					interface.push( name_push, locp, cpp_push[i][j][k] );
				}
            }
        }

        int sent = interface.commit( n );

		// push data to the other solver
		for ( int i = 0; i < Nx; ++i ) {
            for ( int j = 0; j < Ny; ++j ) {
				for ( int k = 0; k < Nz; ++k ) {
					point3d locf( pf[i][j][k][0], pf[i][j][k][1], pf[i][j][k][2] );
					python_fetch[i][j][k] = interface.fetch( name_fetch, locf, 
						n, 
						s1, 
						s2 );
				}
            }
        }

		for ( int i = 0; i < Nx; ++i ) {
			for ( int j = 0; j < Ny; ++j ) {
				for ( int k = 0; k < Nz; ++k ) {
					printf( "{PUSHER_FETCHER_0} python_fetch[%d][%d][%d]: %lf\n", i, j, k, python_fetch[i][j][k] );
				}
			}
		}

	}

    return 0;
}