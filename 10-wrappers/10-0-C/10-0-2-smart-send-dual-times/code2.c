/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
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
 * @file code2.c
 * @author W. Liu
 * @date 24 March 2023
 * @brief C demo to show smart send functions with dual time types
 */

// Standard C includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Include general non-dimensional MUI functions
#include "mui_c_wrapper_general.h"
// Include 3D MUI functions
#include "mui_c_wrapper_3d.h"

int main(int argc, char **argv) {

	// MUI/MPI variables
	double zero = 0.0;
	double tolerance = 8e-1;
	int reset_log = 1;
	int upper_forget = 5;
	int synchronised = 1;
	int mui_ranks, mui_size;
	MPI_Comm MUI_COMM_WORLD;

	// Local parameters
	const int Nx = 2; // number of grid points in x axis
	const int Ny = 2; // number of grid points in y axis
	const int Nz = 2; // number of grid points in z axis
	const int steps = 2; // total time steps
	const int itersteps = 2; // total iteration steps
	const double rSearch = 1.0; // search radius
	char *header = "mpi://";
	char *appname = "PUSHER_FETCHER_2";
	char *interface = "/interface";
	char *footer = "\0";
	char *name_fetch = "code1Value";
	char *name_push = "code2Value";

	int Nt = Nx * Ny * Nz;
	int i, j, k, n, iter;
	double x, y, z;

	double pp[Nx][Ny][Nz][3];
	double pf[Nx][Ny][Nz][3];
	double value_push[Nx][Ny][Nz];
	double value_fetch[Nx][Ny][Nz];

	// Call mui_mpi_split_by_app() function to obtain MPI communicator
	MUI_COMM_WORLD=mui_mpi_split_by_app();

	// Create URI
	char *uri3d = (char*) malloc(strlen(header) + strlen(appname) + strlen(interface) + strlen(footer));
	strcpy(uri3d, header);
	strcat(uri3d, appname);
	strcat(uri3d, interface);
	strcat(uri3d, footer);
	printf("{%s}: URI: %s\n", appname, uri3d);

	// Create MUI interface
	mui_uniface_3d *uniface3d = mui_create_uniface_3d((const char*) uri3d);

    MPI_Comm_rank(MUI_COMM_WORLD, &mui_ranks);
    MPI_Comm_size(MUI_COMM_WORLD, &mui_size);
    printf("{%s}: COMM_SIZE: %d COMM_RANK: %d\n", appname, mui_size, mui_ranks);

	// Create spatial and temporal samplers for fetch operation
    mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler3d = mui_create_sampler_pseudo_nearest_neighbor_3d(rSearch);
	mui_temporal_sampler_exact_3d *temporal_sampler3d = mui_create_temporal_sampler_exact_3d(tolerance);

	// Define bounding box
	double local_x0 = 4.0; // local push origin box start
	double local_y0 = 0.0;
	double local_z0 = 0.0 + ((2.0/((double)mui_size))*((double)mui_ranks));
	double local_x1 = 5.0; // local push origin box end
	double local_y1 = 1.0;
	double local_z1 = 0.0 + ((2.0/((double)mui_size))*((double)mui_ranks+1));
	double local_x2 = 6.0; // local fetch origin box start
	double local_y2 = 0.0;
	double local_z2 = 0.0 + ((2.0/((double)mui_size))*((double)mui_ranks));
	double local_x3 = 7.0; // local fetch origin box end
	double local_y3 = 1.0;
	double local_z3 = 0.0 + ((2.0/((double)mui_size))*((double)mui_ranks+1));

	// Set push points & value
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				double x = local_x0 + (((local_x1 - local_x0)/((double)Nx))*(double)i);
				double y = local_y0 + (((local_y1 - local_y0)/((double)Ny))*(double)j);
				double z = local_z0 + (((local_z1 - local_z0)/((double)Nz))*(double)k);
				pp[i][j][k][0] = x;
				pp[i][j][k][1] = y;
				pp[i][j][k][2] = z;
				value_push[i][j][k] = 97.7777;
			}
		}
	}

	// Set fetch points & value
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				double x = local_x2 + (((local_x3 - local_x2)/((double)Nx))*(double)i);
				double y = local_y2 + (((local_y3 - local_y2)/((double)Ny))*(double)j);
				double z = local_z2 + (((local_z3 - local_z2)/((double)Nz))*(double)k);
				pf[i][j][k][0] = x;
				pf[i][j][k][1] = y;
				pf[i][j][k][2] = z;
				value_fetch[i][j][k] = 11.11;
			}
		}
	}

	// MUI annouce send/rcv span
	mui_announce_send_span_3d_box(uniface3d, local_x0, local_y0, local_z0, local_x1, local_y1, local_z1, zero, steps, synchronised);
	mui_announce_recv_span_3d_box(uniface3d, local_x2, local_y2, local_z2, local_x3, local_y3, local_z3, zero, steps, synchronised);

	for (int n = 1; n <= steps; ++n) {
		for (int iter = 1; iter <= itersteps; ++iter) {

		    printf("{%s}: %d Step %d Sub-iteration\n", appname, n, iter);

			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					for (int k = 0; k < Nz; k++) {
						// Define point for push
						mui_point_3d push_point3d = {pp[i][j][k][0], pp[i][j][k][1], pp[i][j][k][2]};
						// MUI push points
						mui_push_3d(uniface3d, name_push, push_point3d, value_push[i][j][k]);
					}
				}
			}

			// MUI commit
			mui_commit_3d_pair(uniface3d, n, iter);

			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					for (int k = 0; k < Nz; k++) {
						// Define point for fetch
						mui_point_3d fetch_point3d = {pf[i][j][k][0], pf[i][j][k][1], pf[i][j][k][2]};

						// MUI fetch points
						value_fetch[i][j][k] = mui_fetch_pseudo_nearest_neighbor_exact_3d_pair(uniface3d, name_fetch, fetch_point3d, n, iter, spatial_sampler3d,
								temporal_sampler3d);
					}
				}
			}

			if ((n-upper_forget) > 0) {
				mui_forget_upper_3d(uniface3d,(n-upper_forget),reset_log);
			}

			// Print fetched values
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					for (int k = 0; k < Nz; k++) {
					    printf("{%s}: %f\n", appname, value_fetch[i][j][k]);
					}
				}
			}
		}
	}

	// Destroy created 3D MUI objects
	mui_destroy_sampler_pseudo_nearest_neighbor_3d(spatial_sampler3d);
	mui_destroy_temporal_sampler_exact_3d(temporal_sampler3d);
	// Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
	mui_destroy_uniface_3d(uniface3d);

	return 0;
}
