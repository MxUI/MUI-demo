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

/**
 * @file 3D-pseudo-diffusion-coarse.cpp
 * @author W. Liu
 * @date 01 November 2019
 * @brief Coarse (middle) domain of the 3D pseudo diffusion case
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

int main(int argc, char **argv) {

  /// Declare MPI common world with the scope of MUI
  MPI_Comm world = mui::mpi_split_by_app();

  /// Declare MPI ranks and rank size
  int rank, size;
  MPI_Comm_rank(world, &rank);
  MPI_Comm_size(world, &size);

  if (size > 2) {
	  std::cout << "MPI Size larger than 2 does not supported yet." << std::endl;
	          exit(EXIT_FAILURE);
  }

  /// Create rbf matrix folder
  std::string makedirMString = "rbfCoarseMatrix" + std::to_string(rank);
  mkdir(makedirMString.c_str(), 0777);

  /// Define the name of MUI domain
  std::string domain = "coarseDomain";

  /// Define the name of MUI interfaces
  std::vector<std::string> interfaces;
  interfaces.emplace_back("interface2D01");
  interfaces.emplace_back("interface2D02");

  /// Declare MUI objects using MUI configure file
  auto ifs = mui::create_uniface < mui::demo7_config > (domain, interfaces);

  /// Define the name of push/fetch values
  const char *name_pushA = "coarseFieldA";
  const char *name_pushB = "coarseFieldB";
  const char *name_pushC = "coarseFieldC";
  const char *name_fetchA = "fineFieldA";
  const char *name_fetchB = "fineFieldB";
  const char *name_fetchC = "fineFieldC";

  /// Define the forget steps of MUI to reduce the memory
  int forgetSteps = 5;

  /// Define parameters of the RBF sampler
  /// Define the search radius of the RBF sampler
  double rSampler   = 0.8;
  int basisFunc     = 1;
  bool conservative = true;
  double cutoff     = 1e-9;
  bool smoothFunc   = false;
  bool writeMatrix  = true;
  double cgSolveTol	= 1e-6;
  int cgMaxIter     = 500;
  int preconditioner= 0;
  int pouSize    	= 50;
  std::string fileAddress(makedirMString);

  /// Setup diffusion rate
  constexpr static double dr = 0.5;

  /// Setup time steps
  constexpr static int steps = 200;

  /// Setup output interval
  constexpr static int outputInterval = 1;

  /// Domain discretization
  constexpr static int Ntx = 9;          /// total number of grid points in x axis
  constexpr static int Nty = 9;          /// total number of grid points in y axis
  constexpr static int Ntz = 9;          /// total number of grid points in z axis

  int Nx, Ny, Nz;
  Nx = Ntx;                              /// number of grid points in x axis per MPI rank
  Ny = Nty;                              /// number of grid points in y axis per MPI rank

  if (rank < (Ntz % size)) {             /// number of grid points in z axis per MPI rank
	  Nz = Ntz/size + 1;
  } else {
	  Nz = Ntz/size;
  }

  int Nt = Nx * Ny * Nz; /// total number of points per MPI rank

  /// Geometry info
  double lx = 1.0; /// length (x-axis direction) of the geometry
  double ly = 1.0; /// length (y-axis direction) of the geometry
  double lzpr;
  if(size == 2) {
	  if((Ntz%2)==0) {
		  double lpz = 1.0/(static_cast<double>(Ntz)-1.0);
		  double lpz_half = lpz/2.0;
		  lzpr = (1.0 / static_cast<double>(size))-lpz_half; /// length (z-axis direction) per MPI rank of the geometry per MPI rank
	  } else {
		  if(rank == 0){
			  lzpr = 1.0 / static_cast<double>(size); /// length (z-axis direction) per MPI rank of the geometry per MPI rank
		  } else {
			  double lpz = 1.0/(static_cast<double>(Ntz)-1.0);
			  lzpr = (1.0 / static_cast<double>(size))-lpz; /// length (z-axis direction) per MPI rank of the geometry per MPI rank
		  }

	  }
  } else if (size == 1) {
	  lzpr = 1.0;
  }
  double x0 = 1.0; /// origin coordinate (x-axis direction) of the geometry
  double y0 = 0.0; /// origin coordinate (y-axis direction) of the geometry
  double z0pr;
  if(rank == 0) {
	  z0pr = 0.0;
  } else if(rank == 1){
	  if((Ntz%2)==0) {
		  double lpz = 1.0/(static_cast<double>(Ntz)-1.0);
		  double lpz_half = lpz/2.0;
		  z0pr = (1.0 / static_cast<double>(size))+lpz_half; /// origin coordinate (z-axis direction) per MPI rank of the geometry
	  } else {
		  double lpz = 1.0/(static_cast<double>(Ntz)-1.0);
		  z0pr = (1.0 / static_cast<double>(size))+lpz; /// origin coordinate (z-axis direction) per MPI rank of the geometry
	  }
  }

  /// Declare points
  double points[Nx][Ny][Nz][3];

  /// Store point coordinates
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        points[i][j][k][0] = x0 + (lx / (Nx - 1)) * i;
        points[i][j][k][1] = y0 + (ly / (Ny - 1)) * j;
        points[i][j][k][2] = z0pr + (lzpr / (Nz - 1)) * k;
      }
    }
  }

  /// Generate initial pseudo scalar field
  double scalar_field_A[Nx][Ny][Nz];
  double scalar_field_B[Nx][Ny][Nz];
  double scalar_field_C[Nx][Ny][Nz];
  double tolerance = (lx / (Nx - 1)) * 0.5;

  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        scalar_field_A[i][j][k] = 0.0;
        scalar_field_B[i][j][k] = 0.0;
        scalar_field_C[i][j][k] = 0.0;
      }
    }
  }

  /// Declare std::vector to store mui::point2d
  std::vector < mui::point2d > point2dvec;

  /// Store mui::point2d that located in the fetch interface
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        if ((points[i][j][k][0] - x0) <= tolerance) {
          mui::point2d ptf(points[i][j][k][1], points[i][j][k][2]);
          point2dvec.push_back(ptf);
        }
      }
    }
  }

  /// Define and announce MUI send/receive span
  mui::geometry::box<mui::demo7_config> send_region( { y0, z0pr },
      { (y0 + ly), (z0pr + lzpr) });
  ifs[1]->announce_send_span(0, steps, send_region);

  /// Define spatial and temporal samplers
  mui::sampler_rbf<mui::demo7_config> spatial_sampler(rSampler, point2dvec,
       basisFunc, conservative, smoothFunc, writeMatrix, fileAddress,
       cutoff, cgSolveTol, cgMaxIter, pouSize, preconditioner, world);
  mui::temporal_sampler_exact<mui::demo7_config> temporal_sampler;

  /// Commit ZERO step of MUI
  ifs[1]->commit(0);

  /// Output the initial pseudo scalar field
  std::ofstream outputFile;
  std::string filename = "coupling_results" + std::to_string(rank) + "/scalar_field_coarse_0.csv";
  outputFile.open(filename);
  outputFile << "\"X\",\"Y\",\"Z\",\"scalar_field_A\",\"scalar_field_B\",\"scalar_field_C\"\n";
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        outputFile << points[i][j][k][0] << "," << points[i][j][k][1] << ","
            << points[i][j][k][2] << "," << scalar_field_A[i][j][k] << "," << scalar_field_B[i][j][k] << "," << scalar_field_C[i][j][k] << ", \n";
      }
    }
  }
  outputFile.close();

  /// Declare integrations of the pseudo scalar field at boundaries
  double intFaceLD2; /// Integration of the pseudo scalar field at the left boundary of Domain 2
  double intFaceRD2; /// Integration of the pseudo scalar field at the right boundary of Domain 2

  /// Define output files for boundary integrations
  std::ofstream outputIntegration;
  std::string outputIntegrationName = "coupling_results" + std::to_string(rank) + "/faceIntegrationD2.txt";
  outputIntegration.open(outputIntegrationName);
  outputIntegration << "\"t\",\"intFaceLD2\",\"intFaceRD2\"\n";
  outputIntegration.close();

  /// Begin time loops
  for (int t = 1; t <= steps; ++t) {

    printf("\n");
    printf("{Coarse Domain} %d Step \n", t);

    /// Reset boundary integrations
    intFaceLD2 = 0.0;
    intFaceRD2 = 0.0;

    /// Loop over points of Domain 2
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

          /// Loop over left boundary points of Domain 2
          if ((points[i][j][k][0] - x0) <= tolerance) {

            mui::point2d locf(points[i][j][k][1], points[i][j][k][2]);

            /// Fetch data from the other solver
            scalar_field_A[i][j][k] = ifs[0]->fetch(name_fetchA, locf, t,
                spatial_sampler, temporal_sampler);
            scalar_field_B[i][j][k] = ifs[0]->fetch(name_fetchB, locf, t,
                spatial_sampler, temporal_sampler);
            scalar_field_C[i][j][k] = ifs[0]->fetch(name_fetchC, locf, t,
                spatial_sampler, temporal_sampler);

            /// Calculate the integration of left boundary points of Domain 2
            intFaceLD2 += scalar_field_A[i][j][k];

          } else { /// Loop over 'internal' points of Domain 2

            /// Calculate the diffusion of pseudo scalar field of Domain 2
            scalar_field_A[i][j][k] += dr
                * (scalar_field_A[(i - 1)][j][k] - scalar_field_A[i][j][k]);
            scalar_field_B[i][j][k] += dr
                * (scalar_field_B[(i - 1)][j][k] - scalar_field_B[i][j][k]);
            scalar_field_C[i][j][k] += dr
                * (scalar_field_C[(i - 1)][j][k] - scalar_field_C[i][j][k]);

            /// Loop over right boundary points of Domain 2
            if (((x0 + lx) - points[i][j][k][0]) <= (tolerance)) {

              mui::point2d locp(points[i][j][k][1], points[i][j][k][2]);

              /// push data to the other solver
              ifs[1]->push(name_pushA, locp, scalar_field_A[i][j][k]);
              ifs[1]->push(name_pushB, locp, scalar_field_B[i][j][k]);
              ifs[1]->push(name_pushC, locp, scalar_field_C[i][j][k]);

              /// Calculate the integration of right boundary points of Domain 2
              intFaceRD2 += scalar_field_A[i][j][k];
            }
          }
        }
      }
    }

    /// Commit 't' step of MUI
    int sent = ifs[1]->commit(t);

    /// Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
    if ((t - forgetSteps) > 0) {
      ifs[1]->forget(t - forgetSteps);
    }

    /// Output the pseudo scalar field and the boundary integrations
    if ((t % outputInterval) == 0) {

      std::ofstream outputFile;
      std::string filename = "coupling_results" + std::to_string(rank) + "/scalar_field_coarse_"
          + std::to_string(t) + ".csv";
      outputFile.open(filename);
      outputFile << "\"X\",\"Y\",\"Z\",\"scalar_field_A\",\"scalar_field_B\",\"scalar_field_C\"\n";
      for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
          for (int i = 0; i < Nx; ++i) {
            outputFile << points[i][j][k][0] << "," << points[i][j][k][1] << ","
                << points[i][j][k][2] << "," << scalar_field_A[i][j][k] << "," << scalar_field_B[i][j][k] << "," << scalar_field_C[i][j][k] << ", \n";
          }
        }
      }
      outputFile.close();

      outputIntegration.open(outputIntegrationName, std::ofstream::app);
      outputIntegration << std::to_string(t) << "," << intFaceLD2 << ","
          << intFaceRD2 << ", \n";
      outputIntegration.close();

    }
  }

  return 0;
}
