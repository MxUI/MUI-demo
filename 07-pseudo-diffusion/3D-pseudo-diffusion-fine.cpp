/*****************************************************************************	
* Multiscale Universal Interface Code Coupling Library Demo 7                *	
*                                                                            *	
* Copyright (C) 2019 W. Liu                                                  *	
*                                                                            *	
* This software is jointly licensed under the Apache License, Version 2.0    *	
* and the GNU General Public License version 3, you may use it according     *	
* to either.                                                                 *	
*                                                                            *	
* ** Apache License, version 2.0 **                                          *	
*                                                                            *	
* Licensed under the Apache License, Version 2.0 (the "License");            *	
* you may not use this file except in compliance with the License.           *	
* You may obtain a copy of the License at                                    *	
*                                                                            *	
* http://www.apache.org/licenses/LICENSE-2.0                                 *	
*                                                                            *	
* Unless required by applicable law or agreed to in writing, software        *	
* distributed under the License is distributed on an "AS IS" BASIS,          *	
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *	
* See the License for the specific language governing permissions and        *	
* limitations under the License.                                             *	
*                                                                            *	
* ** GNU General Public License, version 3 **                                *	
*                                                                            *	
* This program is free software: you can redistribute it and/or modify       *	
* it under the terms of the GNU General Public License as published by       *	
* the Free Software Foundation, either version 3 of the License, or          *	
* (at your option) any later version.                                        *	
*                                                                            *	
* This program is distributed in the hope that it will be useful,            *	
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *	
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *	
* GNU General Public License for more details.                               *	
*                                                                            *	
* You should have received a copy of the GNU General Public License          *	
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *	
******************************************************************************/	

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

  /// Declare MPI common world with the scope of MUI
  MPI_Comm  world = mui::mpi_split_by_app();

  /// Declare MPI ranks and rank size
  int rank, size;
  MPI_Comm_rank( world, &rank );
  MPI_Comm_size( world, &size );

  /// Create results folder
  std::string makedirCRString = "coupling_results" + std::to_string(rank);
  mkdir(makedirCRString.c_str(), 0777);

  /// Create rbf matrix folder
  std::string makedirMString = "rbfFineMatrix" + std::to_string(rank);
  mkdir(makedirMString.c_str(), 0777);

  /// Define the name of MUI domain
  std::string domain = "fineDomain";

  /// Define the name of MUI interfaces
  std::vector<std::string> interfaces;
  interfaces.emplace_back( "interface2D01" );
  interfaces.emplace_back( "interface2D02" );

  /// Declare MUI objects using MUI configure file
  auto ifs = mui::create_uniface<mui::demo7_config>( domain, interfaces );

  /// Define the name of push/fetch values
  const char* name_pushA   = "fineFieldA";
  const char* name_pushB   = "fineFieldB";
  const char* name_pushC   = "fineFieldC";
  const char* name_fetchA  = "coarseFieldA";
  const char* name_fetchB  = "coarseFieldB";
  const char* name_fetchC  = "coarseFieldC";

  /// Define the forget steps of MUI to reduce the memory
  int forgetSteps = 5;

  /// Define parameters of the RBF sampler
  /// Define the search radius of the RBF sampler
  /// The search radius should not set to a very large value so that to ensure a good convergence
  double rSampler   = 1.0;
  int basisFunc     = 1;
  bool conservative = true;
  double cutoff		= 1e-9;
  bool smoothFunc   = false;
  bool readMatrix  = false;
  bool writeMatrix  = true;
  double cgSolveTol	= 1e-6;
  int cgMaxIter     = 500;
  int preconditioner= 1;
  int pouSize    	= 50;
  std::string fileAddress(makedirMString);

	/// Setup diffusion rate
  constexpr static double dr              = 0.5;

	/// Setup time steps
  constexpr static int    steps           = 200;

  /// Setup output interval
  constexpr static int    outputInterval  = 1;

  /// Geometry info
  double  lx0  = 1.0; /// length (x-axis direction) of the left part geometry
  double  ly0  = 1.0; /// length (y-axis direction) of the left part geometry
  double  lz0pr= 1.0 / static_cast<double>(size); /// length (z-axis direction) per MPI rank of the left part geometry
  double  lz0  = 1.0; /// length (z-axis direction) of the left part geometry
  double  x0   = 0.0;  /// origin coordinate (x-axis direction) of the left part geometry
  double  y0   = 0.0;  /// origin coordinate (y-axis direction) of the left part geometry
  double  z0pr = 0.0 + (lz0pr * static_cast<double>(rank) );  /// origin coordinate (z-axis direction) per MPI rank of the left part geometry
  double  z0   = 0.0;  /// origin coordinate (z-axis direction) of the left part geometry

  double  lx1  = 1.0; /// length (x-axis direction) of the right part geometry
  double  ly1  = 1.0; /// length (y-axis direction) of the right part geometry
  double  lz1pr= 1.0 / static_cast<double>(size); /// length (z-axis direction) per MPI rank of the right part geometry
  double  lz1  = 1.0; /// length (z-axis direction) of the right part geometry
  double  x1   = 2.0;  /// origin coordinate (x-axis direction) of the right part geometry
  double  y1   = 0.0;  /// origin coordinate (y-axis direction) of the right part geometry
  double  z1pr = 0.0 + (lz1pr * static_cast<double>(rank) );  /// origin coordinate (z-axis direction) per MPI rank of the right part geometry
  double  z1   = 0.0;  /// origin coordinate (z-axis direction) of the right part geometry

  /// Domain discretization
  constexpr static int Ntx = 18;         /// number of grid points in x axis per part
  constexpr static int Nty = 18;         /// number of grid points in y axis per part
  constexpr static int Ntz = 18;         /// number of grid points in z axis per part

  int Nx, Ny, Nz;
  Nx = Ntx;                              /// number of grid points in x axis per MPI rank
  Ny = Nty;                              /// number of grid points in y axis per MPI rank

  if (rank < (Ntz % size)) {             /// number of grid points in z axis per MPI rank
	  Nz = Ntz/size + 1;
  } else {
	  Nz = Ntz/size;
  }

  int Nt = 2* Nx * Ny * Nz; /// total number of points per MPI rank

  /// Declare points
  double points0[Nx][Ny][Nz][3], points1[Nx][Ny][Nz][3];

  /// Store point coordinates
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        points0[i][j][k][0] = x0 + (lx0 / (Nx - 1)) * i;
        points0[i][j][k][1] = y0 + (ly0 / (Ny - 1)) * j;

		if(rank==0) {
			points0[i][j][k][2] = z0pr + (lz0pr / (Nz - 1)) * k;
		} else {
			points0[i][j][k][2] = z0pr + (lz0pr / (Nz)) * (k+1);
		}

        points1[i][j][k][0] = x1 + (lx1 / (Nx - 1)) * i;
        points1[i][j][k][1] = y1 + (ly1 / (Ny - 1)) * j;

		if(rank==0) {
			points1[i][j][k][2] = z1pr + (lz1pr / (Nz - 1)) * k;
		} else {
			points1[i][j][k][2] = z1pr + (lz1pr / (Nz)) * (k+1);
		}

      }
    }
  }

  /// Generate initial pseudo scalar fields
  double scalar_fieldA0[Nx][Ny][Nz], scalar_fieldA1[Nx][Ny][Nz];
  double scalar_fieldB0[Nx][Ny][Nz], scalar_fieldB1[Nx][Ny][Nz];
  double scalar_fieldC0[Nx][Ny][Nz], scalar_fieldC1[Nx][Ny][Nz];
  double tolerance = (lx0 / (Nx - 1)) * 0.5;
  double amplitude = 100.0;
  double x_center = x0 + lx0 / 2.0;
  double y_center = y0 + ly0 / 2.0;
  double z_center = z0 + lz0 / 2.0;
  double r_max = sqrt(pow((x0 - x_center), 2) + pow((y0 - y_center), 2)+ pow((z0 - z_center), 2));

  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {

        if (points0[i][j][k][0] <= tolerance) {

          scalar_fieldA0[i][j][k] = amplitude * 2 * (0.5 - (std::max((std::fabs(points0[i][j][k][1]-0.5)),
																	(std::fabs(points0[i][j][k][2]-0.5)))));

          double r_ = sqrt(
              pow((points0[i][j][k][1] - y_center), 2)
                  + pow((points0[i][j][k][2] - z_center), 2));
          scalar_fieldB0[i][j][k] = amplitude * (1 - ((r_max - r_) / r_max));
          scalar_fieldC0[i][j][k] = amplitude * (r_max - r_) / r_max;
        } else {
          scalar_fieldA0[i][j][k] = 0.0;
          scalar_fieldB0[i][j][k] = 0.0;
          scalar_fieldC0[i][j][k] = 0.0;
        }
        scalar_fieldA1[i][j][k] = 0.0;
        scalar_fieldB1[i][j][k] = 0.0;
        scalar_fieldC1[i][j][k] = 0.0;
      }
    }
  }

  /// Declare std::vector to store mui::point2d
  std::vector<mui::point2d> point2dvec;

  /// Store mui::point2d that located in the fetch interface
 	for ( int k = 0; k < Nz; ++k ) {
    for ( int j = 0; j < Ny; ++j ) {
      for ( int i = 0; i < Nx; ++i ) {
        if ((points1[i][j][k][0] - x1) <= tolerance){
          mui::point2d ptf( points1[i][j][k][1], points1[i][j][k][2]);
          point2dvec.push_back(ptf);
        }
      }
    }
  }

  /// Define and announce MUI send/receive span
  mui::geometry::box<mui::demo7_config> send_region( {y0, z0pr}, {(y0+ly0), (z0pr+lz0pr)} );
  mui::geometry::box<mui::demo7_config> recv_region( {y1, z1pr}, {(y1+ly1), (z1pr+lz1pr)} );
  ifs[0]->announce_send_span( 0, steps, send_region );
  ifs[1]->announce_recv_span( 0, steps, recv_region );

	/// Define spatial and temporal samplers
  mui::sampler_rbf<mui::demo7_config> spatial_sampler(rSampler, point2dvec,
      basisFunc, conservative, smoothFunc, writeMatrix, fileAddress,
      cutoff, cgSolveTol, cgMaxIter, pouSize, preconditioner, world);
  mui::temporal_sampler_exact<mui::demo7_config> temporal_sampler;

  /// Commit ZERO step of MUI
  ifs[0]->commit(0);

  /// Output the initial pseudo scalar field
  std::ofstream outputFileLeft;
  std::string filenameL = "coupling_results" + std::to_string(rank) + "/scalar_field_left_fine_0.csv";
  outputFileLeft.open(filenameL);
  outputFileLeft << "\"X\",\"Y\",\"Z\",\"scalar_field_A\",\"scalar_field_B\",\"scalar_field_C\"\n";
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        outputFileLeft << points0[i][j][k][0] << "," << points0[i][j][k][1]
            << "," << points0[i][j][k][2] << "," << scalar_fieldA0[i][j][k] << "," << scalar_fieldB0[i][j][k] << "," << scalar_fieldC0[i][j][k]
            << ", \n";
      }
    }
  }
  outputFileLeft.close();

  std::ofstream outputFileRight;
  std::string filenameR = "coupling_results"  + std::to_string(rank) + "/scalar_field_right_fine_0.csv";
  outputFileRight.open(filenameR);
  outputFileRight << "\"X\",\"Y\",\"Z\",\"scalar_field_A\",\"scalar_field_B\",\"scalar_field_C\"\n";
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        outputFileRight << points1[i][j][k][0] << "," << points1[i][j][k][1]
            << "," << points1[i][j][k][2] << "," << scalar_fieldA1[i][j][k] << "," << scalar_fieldB1[i][j][k] << "," << scalar_fieldC1[i][j][k]
            << ", \n";
      }
    }
  }
  outputFileRight.close();

  /// Declare integrations of the pseudo scalar field at boundaries
  double intFaceLD1; // Integration of the pseudo scalar field at the left boundary of Domain 1
  double intFaceRD1; // Integration of the pseudo scalar field at the right boundary of Domain 1
  double intFaceLD3; // Integration of the pseudo scalar field at the left boundary of Domain 3
  double intFaceRD3; // Integration of the pseudo scalar field at the right boundary of Domain 3

  /// Define output files for boundary integrations
  std::ofstream outputIntegration;
  std::string outputIntegrationName =
      "coupling_results"  + std::to_string(rank) + "/faceIntegrationD1_D3.txt";
  outputIntegration.open(outputIntegrationName);
  outputIntegration
      << "\"t\",\"intFaceLD1\",\"intFaceRD1\",\"intFaceLD3\",\"intFaceRD3\"\n";
  outputIntegration.close();

  /// Begin time loops
  for (int t = 1; t <= steps; ++t) {

    printf("\n");
    printf("{Fine Domain} %d Step \n", t);

    /// Reset boundary integrations
    intFaceLD1 = 0.0;
    intFaceRD1 = 0.0;
    intFaceLD3 = 0.0;
    intFaceRD3 = 0.0;

    /// Loop over points of Domain 1
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

          /// Loop over 'internal' points of Domain 1
          if ((points0[i][j][k][0] - x0) > tolerance) {

            /// Calculate the diffusion of pseudo scalar field of Domain 1
            scalar_fieldA0[i][j][k] += dr
                * (scalar_fieldA0[(i - 1)][j][k] - scalar_fieldA0[i][j][k]);
            scalar_fieldB0[i][j][k] += dr
                * (scalar_fieldB0[(i - 1)][j][k] - scalar_fieldB0[i][j][k]);
            scalar_fieldC0[i][j][k] += dr
                * (scalar_fieldC0[(i - 1)][j][k] - scalar_fieldC0[i][j][k]);

            /// Loop over right boundary points of Domain 1
            if (((x0 + lx0) - points0[i][j][k][0]) <= (tolerance)) {

              mui::point2d locp(points0[i][j][k][1], points0[i][j][k][2]);

              /// push data to the other solver
              ifs[0]->push(name_pushA, locp, scalar_fieldA0[i][j][k]);
              ifs[0]->push(name_pushB, locp, scalar_fieldB0[i][j][k]);
              ifs[0]->push(name_pushC, locp, scalar_fieldC0[i][j][k]);

              /// Calculate the integration of right boundary points of Domain 1
              intFaceRD1 += scalar_fieldA0[i][j][k];
            }

          } else { /// Loop over left boundary points of Domain 1

            /// Calculate the integration of left boundary points of Domain 1
            intFaceLD1 += scalar_fieldA0[i][j][k];

          }
        }
      }
    }

    /// Commit 't' step of MUI
    int sent = ifs[0]->commit(t);

    /// Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
    if ((t - forgetSteps) > 0) {
      ifs[0]->forget(t - forgetSteps);
    }

    /// Loop over points of Domain 3
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

          /// Loop over left boundary points of Domain 3
          if ((points1[i][j][k][0] - x1) <= tolerance) {

            mui::point2d locf(points1[i][j][k][1], points1[i][j][k][2]);

            /// Fetch data from the other solver
            scalar_fieldA1[i][j][k] = ifs[1]->fetch(name_fetchA, locf, t,
                spatial_sampler, temporal_sampler);
            scalar_fieldB1[i][j][k] = ifs[1]->fetch(name_fetchB, locf, t,
                spatial_sampler, temporal_sampler);
            scalar_fieldC1[i][j][k] = ifs[1]->fetch(name_fetchC, locf, t,
                spatial_sampler, temporal_sampler);

            /// Calculate the integration of left boundary points of Domain 3
            intFaceLD3 += scalar_fieldA1[i][j][k];

          } else { /// Loop over 'internal' points of Domain 3

            /// Calculate the diffusion of pseudo scalar field of Domain 3
            scalar_fieldA1[i][j][k] += dr
                * (scalar_fieldA1[(i - 1)][j][k] - scalar_fieldA1[i][j][k]);
            scalar_fieldB1[i][j][k] += dr
                * (scalar_fieldB1[(i - 1)][j][k] - scalar_fieldB1[i][j][k]);
            scalar_fieldC1[i][j][k] += dr
                * (scalar_fieldC1[(i - 1)][j][k] - scalar_fieldC1[i][j][k]);

            /// Loop over right boundary points of Domain 3
            if (((x1 + lx1) - points1[i][j][k][0]) <= (tolerance)) {

              /// Calculate the integration of right boundary points of Domain 3
              intFaceRD3 += scalar_fieldA1[i][j][k];
            }
          }
        }
      }
    }

    /// Output the pseudo scalar field and the boundary integrations
    if ((t % outputInterval) == 0) {

      std::ofstream outputFileLeft;
      std::string filenameL = "coupling_results" + std::to_string(rank) + "/scalar_field_left_fine_"
          + std::to_string(t) + ".csv";
      outputFileLeft.open(filenameL);
      outputFileLeft << "\"X\",\"Y\",\"Z\",\"scalar_field_A\",\"scalar_field_B\",\"scalar_field_C\"\n";

      for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
          for (int i = 0; i < Nx; ++i) {
            outputFileLeft << points0[i][j][k][0] << "," << points0[i][j][k][1]
                << "," << points0[i][j][k][2] << "," << scalar_fieldA0[i][j][k] << "," << scalar_fieldB0[i][j][k] << "," << scalar_fieldC0[i][j][k]
                << ", \n";
          }
        }
      }
      outputFileLeft.close();

      std::ofstream outputFileRight;
      std::string filenameR = "coupling_results" + std::to_string(rank) + "/scalar_field_right_fine_"
          + std::to_string(t) + ".csv";
      outputFileRight.open(filenameR);
      outputFileRight << "\"X\",\"Y\",\"Z\",\"scalar_field_A\",\"scalar_field_B\",\"scalar_field_C\"\n";
      for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
          for (int i = 0; i < Nx; ++i) {
            outputFileRight << points1[i][j][k][0] << "," << points1[i][j][k][1]
                << "," << points1[i][j][k][2] << "," << scalar_fieldA1[i][j][k] << "," << scalar_fieldB1[i][j][k] << "," << scalar_fieldC1[i][j][k]
                << ", \n";
          }
        }
      }
      outputFileRight.close();

      outputIntegration.open(outputIntegrationName, std::ofstream::app);
      outputIntegration << std::to_string(t) << "," << intFaceLD1 << ","
          << intFaceRD1 << "," << intFaceLD3 << "," << intFaceRD3 << ", \n";
      outputIntegration.close();
    }
  }

  return 0;
}
