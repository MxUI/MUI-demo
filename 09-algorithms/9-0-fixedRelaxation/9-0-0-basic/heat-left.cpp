/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 9-0-0            *
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
 * @file heat-left.cpp
 * @author W. Liu
 * @date 16 July 2019
 * @brief Left domain on basic Fixed Relaxation to demonstrate MUI coupling algorithm.
 */

#include "mui.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>

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

    const static int N = 7;
    double u1[N], u2[N];

    /// Initialise values from file
    std::string inoutFilenameL = "Resources/left_FR.csv";
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

        for ( int i = 0; i <  7; i++ ) u1[i] = stod(content[i+1][1]);

    } else {
        std::cerr<<"left_FR.csv missing" << std::endl;
        u1[0] = 1.;
        for ( int i = 1; i <  7; i++ ) u1[i] = 0.;
    }

    inFile.close();

    uniface1d interface( "mpi://left/ifs" );

    MPI_Comm  world = mui::mpi_split_by_app();
    MPI_Comm*  Cppworld = &world;
    int rankLocal = MPI::COMM_WORLD.Get_rank();
    int sizeLocal = MPI::COMM_WORLD.Get_size();
    
    int rank, size;
    MPI_Comm_rank( world, &rank );
    MPI_Comm_size( world, &size );

    /// Create folder
    std::string makedirMString = "results_left" + std::to_string(rank);
    mkdir(makedirMString.c_str(), 0777);
    std::string fileAddress(makedirMString);

    double        k = 0.515, H = 1;
    double *      u = u1, *v = u2;

    std::vector<std::pair<mui::point1d, double>> ptsVluInit;

    for ( int i = 1; i <  7; i++ ) {
        mui::point1d pt(i);
        ptsVluInit.push_back(std::make_pair(pt,u1[i]));
    }

    // fetch data from the other solver
    sampler_pseudo_nearest_neighbor1d<double> s1(0.1);
    temporal_sampler_exact1d  s2;
    algo_fixed_relaxation1d fr(0.01,ptsVluInit);

     // Print off a hello world message
    printf("Hello world from Left rank %d out of %d MUI processors\n",
           rank, size);
           
     // Print off a hello world message
    printf("Hello world from Left rank %d out of %d local processors\n",
           rankLocal, sizeLocal);

    /// Output
    std::ofstream outputFileLeft;
    std::string filenameL = "results_left" + std::to_string(rank) + "/solution-left_FR_0.csv";
    outputFileLeft.open(filenameL);
    outputFileLeft << "\"X\",\"u\"\n";
    for ( int i = 0; i <  7; i++ ) outputFileLeft << i * H << "," << u[i] << ", \n";
    outputFileLeft.close();

    for ( int iter = 1; iter <= 1000; ++iter ) {
        printf( "Left grid itertion %d\n", iter );

            // push data to the other solver
            interface.push( "u", 4, u[4]);
            interface.commit( iter );

            u[6] = interface.fetch( "u0", 6 * H, iter, s1, s2, fr );


            // calculate 'interior' points
            for ( int i = 1; i <  6; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
            // calculate 'boundary' points
            v[0]     = 1.0;

            v[N - 1] = u[N - 1]; 

        // I/O
        std::swap( u, v );
        /// Output
        std::ofstream outputFileLeft;
        std::string filenameL = "results_left" + std::to_string(rank) + "/solution-left_FR_"
          + std::to_string(iter) + ".csv";
        outputFileLeft.open(filenameL);
        outputFileLeft << "\"X\",\"u\"\n";
        for ( int i = 0; i <  7; i++ ) outputFileLeft << i * H << "," << u[i] << ", \n";
        outputFileLeft.close();

    }

    return 0;
}
