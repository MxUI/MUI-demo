/**
 *
 * @file 3D_pseudo_diffusion_coarse.cpp
 * @author W.L.
 * @date 09 Oct 2019
 * @brief See README.md
 *
 */

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

/// Include MUI header file and configure file 
#include "mui.h"
#include "demo7_config.h"

int main(int argc, char ** argv) {

    /// Create rbf matrix folder
    mkdir("rbfCoarseMatrix", 0777);

    /// Declare MPI common world with the scope of MUI
	MPI_Comm  world = mui::mpi_split_by_app();

    /// Declare MPI ranks and rank size
    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    /// Define the name of MUI domain
    std::string domain = "coarseDomain";

    /// Define the name of MUI interfaces 
    std::vector<std::string> interfaces;
    interfaces.emplace_back( "interface2D01" );
    interfaces.emplace_back( "interface2D02" );

    /// Declare MUI objects using MUI configure file
    auto ifs = mui::create_uniface<mui::demo7_config>( domain, interfaces );

    /// Define the name of push/fetch values
	const char* name_push   = "coarseField";			
	const char* name_fetch  = "fineField";

    /// Define the forget steps of MUI to reduce the memory
    int         forgetSteps = 5;

    /// Define parameters of the RBF sampler
    /// Define the search radius of the RBF sampler
    /// The search radius should not set to a very large value so that to ensure a good convergence
    double      rSampler    = 	0.4;
	bool 		conservative=	false; 
	double 		cutoff		=	1e-9;
	bool 		polynomial	=	true;
	bool 		readMatrix	=	false;
    std::string fileAddress("rbfCoarseMatrix");

	/// Setup diffusion rate
    constexpr static double dr              = 0.5;

	/// Setup time steps
    constexpr static int    steps           = 200;

    /// Setup output interval
    constexpr static int    outputInterval  = 1;

    /// Geometry info
    double  x0   = 1.0; /// origin coordinate (x-axis direction) of the geometry
    double  y0   = 0.0; /// origin coordinate (y-axis direction) of the geometry
    double  z0   = 0.0; /// origin coordinate (z-axis direction) of the geometry
    double  lx   = 1.0; /// length (x-axis direction) of the geometry
    double  ly   = 1.0; /// length (y-axis direction) of the geometry
    double  lz   = 1.0; /// length (z-axis direction) of the geometry

    /// Domain discretization
    constexpr static int    Nx      = 9;            /// number of grid points in x axis
    constexpr static int    Ny      = 9;            /// number of grid points in y axis
    constexpr static int    Nz      = 9;            /// number of grid points in z axis
    constexpr static int    Nt      = Nx * Ny * Nz; /// total number of points

    /// Declare points
    double points[Nx][Ny][Nz][3];

    /// Store point coordinates
    for ( int k = 0; k < Nz; ++k ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int i = 0; i < Nx; ++i ) {
                points[i][j][k][0] = x0 + (lx/(Nx-1)) * i;
                points[i][j][k][1] = y0 + (ly/(Ny-1)) * j;
                points[i][j][k][2] = z0 + (lz/(Nz-1)) * k;
            }
        }
    }

    /// Generate initial pseudo scalar field
    double scalar_field[Nx][Ny][Nz];
    double tolerance = (lx/(Nx-1))*0.5;

 	for ( int k = 0; k < Nz; ++k ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int i = 0; i < Nx; ++i ) {
                scalar_field[i][j][k] = 0.0;
            }
        }
    }

    /// Declare std::vector to store mui::point2d
    std::vector<mui::point2d> point2dvec;

    /// Store mui::point2d that located in the fetch interface
 	for ( int k = 0; k < Nz; ++k ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int i = 0; i < Nx; ++i ) {
                if ((points[i][j][k][0] - x0) <= tolerance){
                    mui::point2d ptf( points[i][j][k][1], points[i][j][k][2]);
                    point2dvec.push_back(ptf);
                }
            }
        }
    }

    /// Define and announce MUI send/receive span
    mui::geometry::box<mui::demo7_config> send_region( {y0, z0}, {(y0+ly), (z0+lz)} );
    mui::geometry::box<mui::demo7_config> recv_region( {y0, z0}, {(y0+ly), (z0+lz)} );
    ifs[1]->announce_send_span( 0, steps, send_region );
    ifs[0]->announce_recv_span( 0, steps, recv_region );

	/// Define spatial and temporal samplers
    mui::sampler_rbf<mui::demo7_config> spatial_sampler(rSampler,point2dvec,conservative,cutoff,polynomial,fileAddress,readMatrix);
    mui::chrono_sampler_exact<mui::demo7_config> chrono_sampler;

	/// Commit ZERO step of MUI
	ifs[1]->commit(0);

    /// Output the initial pseudo scalar field
    std::ofstream outputFile;
    std::string filename = "coupling_results/scalar_field_coarse_0.csv";
    outputFile.open (filename);
    outputFile << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
 	for ( int k = 0; k < Nz; ++k ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int i = 0; i < Nx; ++i ) {
                outputFile << points[i][j][k][0] << "," <<points[i][j][k][1]<< "," << points[i][j][k][2]<< "," << scalar_field[i][j][k] << ", \n";
            }
        }
    }
    outputFile.close();

    /// Declare integrations of the pseudo scalar field at boundaries
    double  intFaceLD2; /// Integration of the pseudo scalar field at the left boundary of Domain 2
    double  intFaceRD2; /// Integration of the pseudo scalar field at the right boundary of Domain 2

    /// Define output files for boundary integrations
    std::ofstream outputIntegration;
    std::string outputIntegrationName = "coupling_results/faceIntegrationD2.txt";
    outputIntegration.open (outputIntegrationName);
    outputIntegration << "\"t\",\"intFaceLD2\",\"intFaceRD2\"\n";
    outputIntegration.close();

	/// Begin time loops
    for ( int t = 1; t <= steps; ++t ) {

		printf("\n");
		printf("{Coarse Domain} %d Step \n", t );

        /// Reset boundary integrations
        intFaceLD2   = 0.0;
        intFaceRD2   = 0.0;

        /// Loop over points of Domain 2
        for ( int k = 0; k < Nz; ++k ) {
            for ( int j = 0; j < Ny; ++j ) {
                for ( int i = 0; i < Nx; ++i ) {

                    /// Loop over left boundary points of Domain 2
                    if ((points[i][j][k][0] - x0) <= tolerance){

                        mui::point2d locf( points[i][j][k][1], points[i][j][k][2]);

                        /// Fetch data from the other solver
                        scalar_field[i][j][k] = ifs[0]->fetch(  name_fetch,
                                                                locf, 
                                                                t, 
                                                                spatial_sampler, 
                                                                chrono_sampler);

                        /// Calculate the integration of left boundary points of Domain 2
                        intFaceLD2 += scalar_field[i][j][k];

                    } else{ /// Loop over 'internal' points of Domain 2

                        /// Calculate the diffusion of pseudo scalar field of Domain 2
                        scalar_field[i][j][k] += dr * (scalar_field[(i-1)][j][k] - scalar_field[i][j][k]);

                        /// Loop over right boundary points of Domain 2
                        if (( (x0+lx) - points[i][j][k][0]) <= (tolerance)){

                            mui::point2d locp( points[i][j][k][1], points[i][j][k][2] );

                            /// push data to the other solver
                            ifs[1]->push( name_push, locp, scalar_field[i][j][k] );

                            /// Calculate the integration of right boundary points of Domain 2
                            intFaceRD2 += scalar_field[i][j][k];
                        }
                    }
                }
            }
        }

        /// Commit 't' step of MUI
        int sent = ifs[1]->commit( t );

        /// Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
        if ( (t - forgetSteps) > 0 ) {
            ifs[1]->forget( t - forgetSteps );
        }

        /// Output the pseudo scalar field and the boundary integrations
        if ((t % outputInterval) == 0){

            std::ofstream outputFile;
            std::string filename = "coupling_results/scalar_field_coarse_" + std::to_string(t) + ".csv";
            outputFile.open (filename);
            outputFile << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
            for ( int k = 0; k < Nz; ++k ) {
                for ( int j = 0; j < Ny; ++j ) {
                    for ( int i = 0; i < Nx; ++i ) {
                        outputFile << points[i][j][k][0] << "," <<points[i][j][k][1]<< "," << points[i][j][k][2]<< "," << scalar_field[i][j][k] << ", \n";
                    }
                }
            }
            outputFile.close();

            outputIntegration.open (outputIntegrationName, std::ofstream::app);
            outputIntegration << std::to_string(t) << "," <<intFaceLD2<< "," << intFaceRD2 << ", \n";
            outputIntegration.close();

        }
	}

    return 0;
}