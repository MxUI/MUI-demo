/**
 *
 * @file 3D_pseudo_diffusion_standalone.cpp
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

int main(int argc, char ** argv) {

    /// Create results folder
    mkdir("standalone_results", 0777);

	/// Setup diffusion rate
    constexpr static double dr              = 0.25;

	/// Setup time steps
    constexpr static int    steps           = 200;

    /// Setup output interval
    constexpr static int    outputInterval  = 20;

    /// Geometry info
    double  x0   = 0.0; /// origin coordinate (x-axis direction) of the geometry
    double  y0   = 0.0; /// origin coordinate (y-axis direction) of the geometry
    double  z0   = 0.0; /// origin coordinate (z-axis direction) of the geometry
    double  lx   = 1.0; /// length (x-axis direction) of the geometry
    double  ly   = 1.0; /// length (y-axis direction) of the geometry
    double  lz   = 1.0; /// length (z-axis direction) of the geometry

    /// Domain discretization
    constexpr static int    Nx      = 11;           /// number of grid points in x axis
    constexpr static int    Ny      = 11;           /// number of grid points in y axis
    constexpr static int    Nz      = 11;           /// number of grid points in z axis
    constexpr static int    Nt      = Nx * Ny * Nz; /// total number of points

    /// Declare points
    double points[Nx][Ny][Nz][3];

    /// Store point coordinates
    for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int k = 0; k < Nz; ++k ) {
                points[i][j][k][0] = x0 + (lx/(Nx-1)) * i;
                points[i][j][k][1] = y0 + (ly/(Ny-1)) * j;
                points[i][j][k][2] = z0 + (lz/(Nz-1)) * k;
            }
        }
    }

    /// Generate initial pseudo scalar field
    double scalar_field[Nx][Ny][Nz];
    double tolerance = (lx/(Nx-1))*0.5;
    double amplitude = 100.0;
    double x_center  = x0 + lx/2.0;
    double y_center  = y0 + ly/2.0;
    double z_center  = z0 + lz/2.0;
    double r_max     = std::sqrt(std::pow((x0-x_center),2)+ std::pow((y0-y_center),2)+ std::pow((z0-z_center),2));

 	for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
			for ( int k = 0; k < Nz; ++k ) {

                 if (points[i][j][k][0] <= tolerance){
                   double r_ = std::sqrt(std::pow((points[i][j][k][1]-y_center),2)+
                                     std::pow((points[i][j][k][2]-z_center),2));

                    scalar_field[i][j][k] = amplitude * (r_max - r_)/r_max;
                } else{
                    scalar_field[i][j][k] = 0.0;
                } 
            }
        }
    }

    /// Output the initial pseudo scalar field
    std::ofstream outputFile;
    std::string filename = "standalone_results/scalar_field_standalone_0.csv";
    outputFile.open (filename);
    outputFile << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
    for ( int i = 0; i < Nx; ++i ) {
        for ( int j = 0; j < Ny; ++j ) {
            for ( int k = 0; k < Nz; ++k ) {
                outputFile << points[i][j][k][0] << "," <<points[i][j][k][1]<< "," << points[i][j][k][2]<< "," << scalar_field[i][j][k] << ", \n";
            }
        }
    }
    outputFile.close();

	/// Begin time loops
    for ( int t = 1; t <= steps; ++t ) {

		printf("\n");
		printf("{Coarse Domain} %d Step \n", t );

        /// Loop over points of Domain
        for ( int i = 0; i < Nx; ++i ) {
            for ( int j = 0; j < Ny; ++j ) {
                for ( int k = 0; k < Nz; ++k ) {

                    /// Loop over 'internal' points of Domain
                    if (points[i][j][k][0] > tolerance){

                        /// Calculate the diffusion of pseudo scalar field of Domain
                        scalar_field[i][j][k] += dr * (scalar_field[(i-1)][j][k] - scalar_field[i][j][k]);
                    }
                }
            }
        }

        /// Output the pseudo scalar field
        if ((t % outputInterval) == 0){

            std::ofstream outputFile;
            std::string filename = "standalone_results/scalar_field_standalone_" + std::to_string(t) + ".csv";
            outputFile.open (filename);
            outputFile << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
            for ( int i = 0; i < Nx; ++i ) {
                for ( int j = 0; j < Ny; ++j ) {
                    for ( int k = 0; k < Nz; ++k ) {
                        outputFile << points[i][j][k][0] << "," <<points[i][j][k][1]<< "," << points[i][j][k][2]<< "," << scalar_field[i][j][k] << ", \n";
                    }
                }
            }
            outputFile.close();
        }
	}

    return 0;
}