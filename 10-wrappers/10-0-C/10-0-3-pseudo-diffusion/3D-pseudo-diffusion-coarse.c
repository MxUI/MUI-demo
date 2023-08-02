/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 10-0-3           *
*                                                                            *
* Copyright (C) 2023 W. Liu                                                  *
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
*****************************************************************************/

/**
 * @file 3D-pseudo-diffusion-coarse.c
 * @author W. Liu
 * @date 04 April 2023
 * @brief Coarse (middle) domain of the 3D pseudo diffusion case
 */

// Standard C includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

// Include general non-dimensional MUI functions
#include "mui_c_wrapper_general.h"
// Include 2D MUI functions
#include "mui_c_wrapper_2d.h"

int main(int argc, char **argv) {

    // MUI/MPI variables
    int mui_ranks, mui_size;
    MPI_Comm MUI_COMM_WORLD;

    // Call mui_mpi_split_by_app() function to obtain MPI communicator
    MUI_COMM_WORLD=mui_mpi_split_by_app();

    /// Declare MPI ranks and rank size
    MPI_Comm_rank(MUI_COMM_WORLD, &mui_ranks);
    MPI_Comm_size(MUI_COMM_WORLD, &mui_size);

    if (mui_size > 2) {
    printf("MPI Size larger than 2 does not supported yet.\n");
    exit(EXIT_FAILURE);
    }

    /// Create rbf matrix folder
    char makedirMString[32];
    sprintf(makedirMString, "rbfCoarseMatrix%d", mui_ranks);
    mkdir(makedirMString, 0777);

    /// Define the name of MUI domain
    char *domain = (char*) malloc(strlen("coarseDomain") + 1);
    strcpy(domain, "coarseDomain");
    strcat(domain, "\0");

    // Get the number of interfaces to create
    int num_interfaces = 2;

    /// Define the name of MUI interfaces
    char** interfaces = (char**) malloc(sizeof(char*) * num_interfaces);
    for(int i=0; i<num_interfaces; i++) {
        interfaces[i] = (char*) malloc(strlen("interface2D0" + 2));
        strcpy(interfaces[i], "interface2D0");
        char tmp[11];
        sprintf(tmp, "%d", (i+1));
        strcat(interfaces[i], tmp);
        strcat(interfaces[i], "\0");
    }
    // Define return structure to hold created unifaces
    mui_uniface_2d** ifs = (mui_uniface_2d**) malloc(sizeof(mui_uniface_2d*) * num_interfaces);

    // Create MUI interfaces using helper function
    ifs = mui_create_uniface_multi_2d( (const char*)domain, (const char**)interfaces, num_interfaces );

    // Free char memory after uniface creation
    for(int i=0; i<num_interfaces; i++) {
        free(interfaces[i]);
    }

    free(domain);

    /// Define the name of push/fetch values
    const char *name_pushA = "coarseFieldA";
    const char *name_fetchA = "fineFieldA";

    /// Define the forget steps of MUI to reduce the memory
    int forgetSteps = 5;

    /// Define parameters of the RBF sampler
    /// Define the search radius of the RBF sampler
    double rSampler   = 0.8;
    int basisFunc     = 1;
    int conservative  = 1;
    double cutoff     = 1e-9;
    int smoothFunc    = 0;
    int writeMatrix   = 1;
    double cgSolveTol = 1e-6;
    int cgMaxIter     = 500;
    int preconditioner= 1;
    int pouSize       = 50;
    int synchronised  = 1;
    int reset_log     = 1;
    char fileAddress[100];
    sprintf(fileAddress, "%s", makedirMString);

    /// Setup diffusion rate
    const double dr = 0.5;

    /// Setup time steps
    const int steps = 200;

    /// Setup output interval
    const int outputInterval = 1;

    /// Domain discretization
    const int Ntx = 9;          /// total number of grid points in x axis
    const int Nty = 9;          /// total number of grid points in y axis
    const int Ntz = 9;          /// total number of grid points in z axis

    int Nx, Ny, Nz;
    Nx = Ntx;                              /// number of grid points in x axis per MPI rank
    Ny = Nty;                              /// number of grid points in y axis per MPI rank

    if (mui_ranks < (Ntz % mui_size)) {    /// number of grid points in z axis per MPI rank
        Nz = Ntz/mui_size + 1;
    } else {
        Nz = Ntz/mui_size;
    }

    int Nt = Nx * Ny * Nz; /// total number

    /// Geometry info
    double lx = 1.0; /// length (x-axis direction) of the geometry
    double ly = 1.0; /// length (y-axis direction) of the geometry
    double lzpr;
    if(mui_size == 2) {
        if((Ntz%2)==0) {
            double lpz = 1.0/((double)(Ntz)-1.0);
            double lpz_half = lpz/2.0;
            lzpr = (1.0 / (double)(mui_size))-lpz_half; /// length (z-axis direction) per MPI rank of the geometry per MPI rank
        } else {
            if(mui_ranks == 0){
                lzpr = 1.0 / (double)(mui_size); /// length (z-axis direction) per MPI rank of the geometry per MPI rank
            } else {
                double lpz = 1.0/((double)(Ntz)-1.0);
                lzpr = (1.0 / (double)(mui_size))-lpz; /// length (z-axis direction) per MPI rank of the geometry per MPI rank
            }
        }
    } else if (mui_size == 1) {
        lzpr = 1.0;
    }
    double x0 = 1.0; /// origin coordinate (x-axis direction) of the geometry
    double y0 = 0.0; /// origin coordinate (y-axis direction) of the geometry
    double z0pr;
    if(mui_ranks == 0) {
        z0pr = 0.0;
    } else if(mui_ranks == 1){
        if((Ntz%2)==0) {
            double lpz = 1.0/((double)(Ntz)-1.0);
            double lpz_half = lpz/2.0;
            z0pr = (1.0 / (double)(mui_size))+lpz_half; /// origin coordinate (z-axis direction) per MPI rank of the geometry
        } else {
            double lpz = 1.0/((double)(Ntz)-1.0);
            z0pr = (1.0 / (double)(mui_size))+lpz; /// origin coordinate (z-axis direction) per MPI rank of the geometry
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
    double tolerance = (lx / (Nx - 1)) * 0.5;

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                scalar_field_A[i][j][k] = 0.0;
            }
        }
    }

    /// Declare array to store mui::point2d
    mui_point_2d point2d_temp[Nt];
    int point_count = 0;

    /// Store mui::point2d that located in the fetch interface
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if ((points[i][j][k][0] - x0) <= tolerance) {
                    point2d_temp[point_count].point_1 = points[i][j][k][1];
                    point2d_temp[point_count].point_2 = points[i][j][k][2];
                    point_count = point_count + 1;
                }
            }
        }
    }

    mui_point_2d point2d[point_count];

    for (int i = 0; i < point_count; ++i) {
        point2d[i].point_1 = point2d_temp[i].point_1;
        point2d[i].point_2 = point2d_temp[i].point_2;
    }

    /// Define spatial and temporal samplers
    mui_sampler_rbf_2d *spatial_sampler2d = mui_create_sampler_rbf_2d(rSampler, point2d, point_count, basisFunc, conservative,
            smoothFunc, writeMatrix, fileAddress, cutoff, cgSolveTol, cgMaxIter, pouSize, preconditioner, MUI_COMM_WORLD);
    mui_temporal_sampler_exact_2d *temporal_sampler2d = mui_create_temporal_sampler_exact_2d(8e-1);

    /// Commit ZERO step of MUI
    mui_commit_2d(ifs[1], 0);

    /// Output the initial pseudo scalar field
    FILE *outputFile;
    char filename[50];
    sprintf(filename, "coupling_results%d/scalar_field_coarse_0.csv", mui_ranks);
    outputFile = fopen(filename, "w");
    fprintf(outputFile, "\"X\",\"Y\",\"Z\",\"scalar_field_A\"\n");
    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                fprintf(outputFile, "%f,%f,%f,%f,\n", points[i][j][k][0], points[i][j][k][1], points[i][j][k][2], scalar_field_A[i][j][k]);
            }
        }
    }
    fclose(outputFile);

    /// Declare integrations of the pseudo scalar field at boundaries
    double intFaceLD2; /// Integration of the pseudo scalar field at the left boundary of Domain 2
    double intFaceRD2; /// Integration of the pseudo scalar field at the right boundary of Domain 2

    /// Define output files for boundary integrations
    FILE *outputIntegration;
    char outputIntegrationName[100];
    sprintf(outputIntegrationName, "coupling_results%d/faceIntegrationD2.txt", mui_ranks);
    outputIntegration = fopen(outputIntegrationName, "w");
    fprintf(outputIntegration, "\"t\",\"intFaceLD2\",\"intFaceRD2\"\n");
    fclose(outputIntegration);

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

                        mui_point_2d locf = {points[i][j][k][1], points[i][j][k][2]};

                        /// Fetch data from the other solver
                        scalar_field_A[i][j][k] = mui_fetch_rbf_exact_2d(ifs[0], name_fetchA,
                            locf, t, spatial_sampler2d, temporal_sampler2d);

                        /// Calculate the integration of left boundary points of Domain 2
                        intFaceLD2 = intFaceLD2 + scalar_field_A[i][j][k];

                    } else { /// Loop over 'internal' points of Domain 2

                        /// Calculate the diffusion of pseudo scalar field of Domain 2
                        scalar_field_A[i][j][k] = scalar_field_A[i][j][k] + dr
                            * (scalar_field_A[(i - 1)][j][k] - scalar_field_A[i][j][k]);

                        /// Loop over right boundary points of Domain 2
                        if (((x0 + lx) - points[i][j][k][0]) <= (tolerance)) {

                            mui_point_2d locp= {points[i][j][k][1], points[i][j][k][2]};

                            /// push data to the other solver
                            mui_push_2d(ifs[1], name_pushA, locp, scalar_field_A[i][j][k]);

                            /// Calculate the integration of right boundary points of Domain 2
                            intFaceRD2 = intFaceRD2 + scalar_field_A[i][j][k];
                        }
                    }
                }
            }
        }

        /// Commit 't' step of MUI
        mui_commit_2d(ifs[1], t);

        /// Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
        if ((t - forgetSteps) > 0) {
            mui_forget_upper_2d(ifs[1],(t - forgetSteps),reset_log);
        }

        /// Output the pseudo scalar field and the boundary integrations
        if ((t % outputInterval) == 0) {
            FILE *outputFile;
            char filename[100];
            sprintf(filename, "coupling_results%d/scalar_field_coarse_%d.csv", mui_ranks, t);
            outputFile = fopen(filename, "w");
            fprintf(outputFile, "\"X\",\"Y\",\"Z\",\"scalar_field_A\"\n");
            for (int k = 0; k < Nz; ++k) {
                for (int j = 0; j < Ny; ++j) {
                    for (int i = 0; i < Nx; ++i) {
                        fprintf(outputFile, "%f,%f,%f,%f,\n",
                                points[i][j][k][0], points[i][j][k][1], points[i][j][k][2],
                                scalar_field_A[i][j][k]);
                    }
                }
            }
            fclose(outputFile);

            FILE *outputIntegration;
            char outputIntegrationName[100];
            sprintf(outputIntegrationName, "coupling_results%d/faceIntegrationD2.txt", mui_ranks);
            outputIntegration = fopen(outputIntegrationName, "a");
            fprintf(outputIntegration, "%d,%f,%f,\n", t, intFaceLD2, intFaceRD2);
            fclose(outputIntegration);
        }
    }

	// Destroy created 3D MUI objects
	mui_destroy_sampler_rbf_2d(spatial_sampler2d);
	mui_destroy_temporal_sampler_exact_2d(temporal_sampler2d);
	// Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
	for(int i=0; i<num_interfaces; i++) {
		mui_destroy_uniface_2d(ifs[i]);
	}

    return 0;
}
