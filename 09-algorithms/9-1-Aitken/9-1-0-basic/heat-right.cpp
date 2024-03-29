/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 9-1-0            *
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
******************************************************************************/

/**
 * @file heat-right.cpp
 * @author W. Liu
 * @date 16 July 2019
 * @brief Right domain on basic Aitken's to demonstrate MUI coupling algorithm.
 */

#include "mui.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

/*                             Left Domain                       Right Domain             
 * Coarse : +-------+-------+-------+-------o=======+=======o-------+-------+-------+-------+
 *          0       1       2       3       4       5       6       7       8       9      10
 * +: grid points
 * o: interface points
 * -: single domain zone
 * =: overlapping zone
 */

// USAGE: mpirun -np 1 ./heat-left : -np 1 ./heat-right

int main( int argc, char ** argv ) {
    using namespace mui;

    const static int N = 11;
    double u1[N], u2[N]; 

    /// Initialise values from file
    std::string inoutFilenameL = "Resources/right_AITKEN.csv";
    std::fstream inFile;
    std::vector<std::vector<std::string>> content;
    std::vector<std::string> row;
    std::string line, word;
    inFile.open(inoutFilenameL, std::ios::in);
    if(inFile.is_open())
    {
        while(getline(inFile, line)) {
            row.clear();
            std::stringstream str(line);
            while (std::getline(str, word, ',')) {
                row.push_back(word);
           }
           content.push_back(row);
        }

        for ( int i = 4; i <  11; i++ ) u1[i] = stod(content[i-3][1]);

    } else {
        std::cerr<<"right_AITKEN.csv missing" << std::endl;
        for ( int i = 4; i <  11; i++ ) u1[i] = 0.0;
    }

    inFile.close();

    uniface1d interface( "mpi://right/ifs" );

    MPI_Comm  world = mui::mpi_split_by_app();
    MPI_Comm*  Cppworld = &world;
    int rankLocal = MPI::COMM_WORLD.Get_rank();
    int sizeLocal = MPI::COMM_WORLD.Get_size();
    
    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    /// Create folder
    std::string makedirMString = "results_right" + std::to_string(rank);
    mkdir(makedirMString.c_str(), 0777);
    std::string fileAddress(makedirMString);
    std::string makedirMIterString = "results_iteration_right" + std::to_string(rank);
    mkdir(makedirMIterString.c_str(), 0777);
    std::string fileAddressIter(makedirMIterString);


    double        k = 0.515, H = 1;
    double *      u = u1, *v = u2;

    std::vector<std::pair<mui::point1d, double>> ptsVluInit;

    for ( int i = 4; i <  11; i++ ) {
        mui::point1d pt(i);
        ptsVluInit.push_back(std::make_pair(pt,u1[i]));
    }

    // fetch data from the other solver
    sampler_pseudo_nearest_neighbor1d<double> s1(0.1);
    temporal_sampler_exact1d  s2;
    algo_aitken1d aitken(0.01,1.0);

     // Print off a hello world message
    printf("Hello world from Right rank %d out of %d MUI processors\n",
           rank, size);
           
     // Print off a hello world message
    printf("Hello world from Right rank %d out of %d local processors\n",
           rankLocal, sizeLocal);

    /// Output
    std::ofstream outputFileRight;
    std::string filenameR = "results_right" + std::to_string(rank) + "/solution-right_AITKEN_0.csv";
    outputFileRight.open(filenameR);
    outputFileRight << "\"X\",\"u\"\n";
    for ( int i = 4; i <  11; i++ ) outputFileRight << i * H << "," << u[i] << ", \n";
    outputFileRight.close();

    std::ofstream outputFileIterRight;
    std::string filenameIterR = "results_iteration_right" + std::to_string(rank) + "/solution-right_AITKEN_0.csv";
    outputFileIterRight.open(filenameIterR);
    outputFileIterRight << "\"X\",\"u\"\n";
    for ( int i = 4; i <  11; i++ ) outputFileIterRight << i * H << "," << u[i] << ", \n";
    outputFileIterRight.close();
    for ( int t = 1; t <= 10; ++t ) {
		for ( int iter = 1; iter <= 100; ++iter ) {
			printf( "Right grid time %d iteration %d\n", t, iter );

				u[4] = interface.fetch( "u", 4 * H,  t, iter, s1, s2 , aitken);
				printf( "Right under relaxation factor at t= %d iter= %d is %f\n", t, iter, aitken.get_under_relaxation_factor(t,iter));
				printf( "Right residual L2 Norm at t= %d iter= %d is %f\n", t, iter, aitken.get_residual_L2_Norm(t,iter));
				// calculate 'interior' points
				for ( int i = 5; i <  11; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
				// calculate 'boundary' points
				v[N - 1] = 0.0;
				v[4]     = u[4    ];

				// push data to the other solver
				interface.push( "u0", 6 * H, u[6] );
				interface.commit(  t, iter );
			// I/O
			std::swap( u, v );

			/// Output
			std::ofstream outputFileIterRight;
			std::string filenameIterR = "results_iteration_right" + std::to_string(rank) + "/solution-right_AITKEN_"
			  + std::to_string(((t-1)*100) + iter) + ".csv";
			outputFileIterRight.open(filenameIterR);
			outputFileIterRight << "\"X\",\"u\"\n";
			for ( int i = 4; i <  11; i++ ) outputFileIterRight << i * H << "," << u[i] << ", \n";
			outputFileIterRight.close();
		}
		/// Output
		std::ofstream outputFileRight;
		std::string filenameR = "results_right" + std::to_string(rank) + "/solution-right_AITKEN_"
		  + std::to_string(t) + ".csv";
		outputFileRight.open(filenameR);
		outputFileRight << "\"X\",\"u\"\n";
		for ( int i = 4; i <  11; i++ ) outputFileRight << i * H << "," << u[i] << ", \n";
		outputFileRight.close();
    }

    return 0;
}
