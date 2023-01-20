/*
 * heat-left.cpp
 *
 *  Created on: Jul 16, 2019
 *      Author: W. Liu
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
    std::string inoutFilenameL = "left_AITKEN.csv";
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
        std::cerr<<"left_AITKEN.csv missing" << std::endl;
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

    /// Create rbf matrix folder
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
    chrono_sampler_exact1d  s2;
    algo_aitken1d aitken(0.01,1.0,world,ptsVluInit,0.023969);

     // Print off a hello world message
    printf("Hello world from Left rank %d out of %d MUI processors\n",
           rank, size);
           
     // Print off a hello world message
    printf("Hello world from Left rank %d out of %d local processors\n",
           rankLocal, sizeLocal);

    /// Output
    std::ofstream outputFileLeft;
    std::string filenameL = "results_left" + std::to_string(rank) + "/solution-left_AITKEN_0.csv";
    outputFileLeft.open(filenameL);
    outputFileLeft << "\"X\",\"u\"\n";
    for ( int i = 0; i <  7; i++ ) outputFileLeft << i * H << "," << u[i] << ", \n";
    outputFileLeft.close();

    for ( int t = 1; t <= 1000; ++t ) {
        printf( "Left grid step %d\n", t );

            // push data to the other solver
            interface.push( "u", 4, u[4]);
            interface.commit( t );

            u[6] = interface.fetch( "u0", 6 * H, t, s1, s2, aitken );

            printf( "Left under relaxation factor at t= %d is %f\n", t, aitken.get_under_relaxation_factor(t));
            printf( "Left residual L2 Norm at t= %d is %f\n", t, aitken.get_residual_L2_Norm(t));

            // calculate 'interior' points
            for ( int i = 1; i <  6; i++ ) v[i] = u[i] + k / ( H * H ) * ( u[i - 1] + u[i + 1] - 2 * u[i] );
            // calculate 'boundary' points
            v[0]     = 1.0;

            v[N - 1] = u[N - 1]; 

        // I/O
        std::swap( u, v );
        /// Output
        std::ofstream outputFileLeft;
        std::string filenameL = "results_left" + std::to_string(rank) + "/solution-left_AITKEN_"
          + std::to_string(t) + ".csv";
        outputFileLeft.open(filenameL);
        outputFileLeft << "\"X\",\"u\"\n";
        for ( int i = 0; i <  7; i++ ) outputFileLeft << i * H << "," << u[i] << ", \n";
        outputFileLeft.close();

    }

    return 0;
}
