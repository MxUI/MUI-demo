/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 10-0-4           *
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
 * @file heat-right.c
 * @author W. Liu
 * @date 25 March 2023
 * @brief Right domain on basic Fixed Relaxation to demonstrate MUI coupling algorithm
 */

// Standard C includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

/*                             Left Domain                       Right Domain
 * Coarse : +-------+-------+-------+-------o=======+=======o-------+-------+-------+-------+
 *          0       1       2       3       4       5       6       7       8       9      10
 * +: grid points
 * o: interface points
 * -: single domain zone
 * =: overlapping zone
 */

// Include general non-dimensional MUI functions
#include "mui_c_wrapper_general.h"
// Include 3D MUI functions
#include "mui_c_wrapper_1d.h"

int main(int argc, char **argv) {

    const static int N = 11;
    double u1[N], u2[N];
    double k = 0.515, H = 1;
    double *u = u1, *v = u2;
	char *header = "mpi://";
	char *appname = "right";
	char *interface = "/ifs";
	char *footer = "\0";
	char *fetch_name = "u";
	char *push_name = "u0";
	const double rSearch = 1.0; // search radius
	double tolerance = 8e-1;
	MPI_Comm MUI_COMM_WORLD;
	int mui_ranks, mui_size;
	mui_point_1d points[N];
	double value_init[N];

    // Open the file in read mode
	FILE* fp = fopen("Resources/right_FR.csv", "r");

    // Check if the file was opened successfully
    if (fp == NULL) {
        printf("right_FR.csv missing.\n");
        return 1;
    }

    // Read each line of the file
    char line[100];
    char* value;
    char* values[10];
    int count;
    char content[100][10][100];
    int i = 0;
    while (fgets(line, sizeof(line), fp) && i < 100) {
        count = 0;
        value = strtok(line, ",");
        while (value != NULL && count < 10) {
            strcpy(content[i][count], value);
            count++;
            value = strtok(NULL, ",");
        }
        i++;
    }

	for ( int i = 4; i <  11; i++ ) u1[i] = atof(content[i-3][1]);

    // Close the file
    fclose(fp);

	// Create URI
	char *uri1d = (char*) malloc(strlen(header) + strlen(appname) + strlen(interface) + strlen(footer));
	strcpy(uri1d, header);
	strcat(uri1d, appname);
	strcat(uri1d, interface);
	strcat(uri1d, footer);
	printf("{%s}: URI: %s\n", appname, uri1d);

	// Create MUI interface
	mui_uniface_1d *uniface1d = mui_create_uniface_1d((const char*) uri1d);

	// Call mui_mpi_split_by_app() function to obtain MPI communicator
	MUI_COMM_WORLD=mui_mpi_split_by_app();

    MPI_Comm_rank(MUI_COMM_WORLD, &mui_ranks);
    MPI_Comm_size(MUI_COMM_WORLD, &mui_size);
    printf("{%s}: COMM_SIZE: %d COMM_RANK: %d\n", appname, mui_size, mui_ranks);

    char makedirMString[32];
    sprintf(makedirMString, "results_right%d", mui_ranks);
    mkdir(makedirMString, 0777);
    char fileAddress[32];
    strcpy(fileAddress, makedirMString);

    int pair_count = 0;
    for ( int i = 4; i <  11; i++ ) {
    	points[i].point_1 = i;
    	value_init[i] = u1[i];
    	pair_count += 1;
    }

	// Create spatial and temporal samplers for fetch operation
    mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler1d = mui_create_sampler_pseudo_nearest_neighbor_1d(10.0);
	mui_temporal_sampler_exact_1d *temporal_sampler1d = mui_create_temporal_sampler_exact_1d(1.0);
	mui_algorithm_fixed_relaxation_1d *algorithm1d = mui_create_algorithm_fixed_relaxation_1d(0.01, MUI_COMM_WORLD, points, value_init, pair_count);

	// Output
	char fileName[100];
	sprintf(fileName, "results_right%d/solution-right_FR_0.csv", mui_ranks);
	FILE* outputFile = fopen(fileName, "w");
	fprintf(outputFile, "\"X\",\"u\"\n");
    for (int i = 4; i < 11; i++) {
        fprintf(outputFile, "%.2f,%.2f,\n", i * H, u[i]);
    }
    fclose(outputFile);

	for (int iter = 1; iter <= 1000; ++iter) {
		printf( "Right grid iteration %d\n", iter );

		// MUI fetch points
		mui_point_1d point_fetch={4 * H};
		u[4] = mui_fetch_pseudo_nearest_neighbor_exact_fixed_relaxation_1d(uniface1d, fetch_name, point_fetch, iter, spatial_sampler1d, temporal_sampler1d, algorithm1d);

		printf( "Right under relaxation factor at iter= %d is %f\n", iter, mui_fixed_relaxation_get_under_relaxation_factor_1d(algorithm1d,iter));
		printf( "Right residual L2 Norm at iter= %d is %f\n", iter, mui_fixed_relaxation_get_residual_L2_Norm_1d(algorithm1d,iter));

        // calculate 'interior' points
        for ( int i = 5; i <  11; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
        // calculate 'boundary' points
        v[N - 1] = 0.0;
        v[4]     = u[4    ];

        // MUI push points
        mui_point_1d point_push;
        point_push.point_1 = 6 * H;
        mui_push_1d(uniface1d, push_name, point_push, u[6]);

		// MUI commit
		mui_commit_1d(uniface1d, iter);

		double* temp = u;
		u = v;
		v = temp;

        // Output
		char fileName[100];
		sprintf(fileName, "results_right%d/solution-right_FR_0%d.csv", mui_ranks,iter);
		FILE* outputFile = fopen(fileName, "w");
		fprintf(outputFile, "\"X\",\"u\"\n");
	    for (int i = 4; i < 11; i++) {
	        fprintf(outputFile, "%.2f,%.2f,\n", i * H, u[i]);
	    }
	    fclose(outputFile);
	}

	// Destroy created 3D MUI objects
	mui_destroy_sampler_pseudo_nearest_neighbor_1d(spatial_sampler1d);
	mui_destroy_temporal_sampler_exact_1d(temporal_sampler1d);
	// Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
	mui_destroy_uniface_1d(uniface1d);

	return 0;
}
