/*
 * heat-right.cpp
 *
 *  Created on: Jul 16, 2019
 *      Author: W. Liu
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

    const static int N = 110;
    double u1[N], u2[N]; 

    /// Initialise values from file
    std::string inoutFilenameL = "right_AITKEN.csv";
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

        for ( int i = 40; i <  110; i+=10 ) u1[i] = stod(content[i*0.1-3][1]);

    } else {
        std::cerr<<"right_AITKEN.csv missing" << std::endl;
        for ( int i = 40; i <  110; i+=10 ) u1[i] = 0.0;
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

    /// Create rbf matrix folder
    std::string makedirMString = "results_right" + std::to_string(rank);
    mkdir(makedirMString.c_str(), 0777);
    std::string fileAddress(makedirMString);

    double        k = 0.515, H = 1;
    double *      u = u1, *v = u2;

    std::vector<std::pair<mui::point1d, double>> ptsVluInit;

    for ( int i = 40; i <  110; i+=10 ) {
        mui::point1d pt(i);
        ptsVluInit.push_back(std::make_pair(pt,u1[i]));
    }

    // fetch data from the other solver
    sampler_pseudo_nearest_neighbor1d<double> s1(30);
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
    for ( int i = 40; i <  110; i+=10 ) outputFileRight << i * H << "," << u[i] << ", \n";
    outputFileRight.close();

    for ( int iter = 1; iter <= 1000; ++iter ) {
        printf( "Right grid iteration %d\n", iter );

            u[40] = interface.fetch( "u", 40 * H, std::numeric_limits<double>::lowest(), iter, s1, s2, aitken );

			if ((iter>=10) && (iter<50)) {
				u[42] = interface.fetch( "u", 42 * H, std::numeric_limits<double>::lowest(), iter, s1, s2, aitken  );
			}

            printf( "Right under relaxation factor at iter= %d is %f\n", iter, aitken.get_under_relaxation_factor(std::numeric_limits<double>::lowest(), iter));
            printf( "Right residual L2 Norm at iter= %d is %f\n", iter, aitken.get_residual_L2_Norm(std::numeric_limits<double>::lowest(), iter));
   
         // calculate 'interior' points
            for ( int i = 50; i <  110; i+=10 ) v[i] = u[i] + k / ( H * H ) * ( u[i - 10] + u[i + 10] - 2 * u[i] );
            // calculate 'boundary' points
            v[N - 10] = 0.0;
            v[40]     = u[40    ];

			if ((iter>=10) && (iter<50)) {
				v[42]     = u[42    ];
			}

            // push data to the other solver
            interface.push( "u0", 60 * H, u[60] );
            interface.commit( std::numeric_limits<double>::lowest(), iter );
        // I/O
        std::swap( u, v );

        /// Output
        std::ofstream outputFileRight;
        std::string filenameR = "results_right" + std::to_string(rank) + "/solution-right_AITKEN_"
          + std::to_string(iter) + ".csv";
        outputFileRight.open(filenameR);
        outputFileRight << "\"X\",\"u\"\n";

		outputFileRight << 40 * H << "," << u[40] << ", \n";

		if ((iter>=10) && (iter<50)) {
			outputFileRight << 42 * H << "," << u[42] << ", \n";
		}

        for ( int i = 50; i <  110; i+=10 ) outputFileRight << i * H << "," << u[i] << ", \n";
        outputFileRight.close();
    }

    return 0;
}
