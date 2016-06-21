// one dimension heat conductivity by SPH
// Author: Xin Bian, 23:53:00, Sunday, June 19, 2016, from Munich.
// Simplification and MUI integration by Yu-Hang Tang

#include <cstdio>
#include <fstream>
#include "../mui/mui.h"

/* Hybrid Lagrangian-Eulerian Simulation of Heat Conduction
 * SPH-Finite Difference, PBC
 * SPH :                         o o o o o o o o o o o o o o o o o o o o o
 *                               0 1 2 3 4 5 6 7 8 9      .....         24
 *                               ^ ^   * *                       * *   ^ ^
 * FDM :    +---+---+---+---+---+---+---+                         +---+---+---+---+---+---+---+
 *          0   1   2   3   4   5   6   7                         8   9  10  11  12  13  14  15
 *                              *   *   ^                         ^   *   *
 * +: grid points
 * o: SPH particles
 * ^: interface points - fetch
 * *: interface points - push
 */

double cubic_spline_gradient( double r, double rc ) {
    const static double sigma  = 2.0 / 3.0; // for one diemension
    double h      = rc / 2.0;
    double coef_g = sigma / h / h;
    double s      = r / h;
    double w_g    = 0;
    if ( s < 1.0 ) w_g = -3.0 * s + 9.0 / 4.0 * s * s;
    else if ( s < 2.0 ) w_g = -3.0 / 4.0 * ( 2 - s ) * ( 2 - s );
    return w_g * coef_g;
}

int main() {

    using namespace mui;
    uniface1d interface( "mpi://sph/ifs" );

    const int Ni = 21, No = 2;//number of inside/outside particles
    const int N = Ni + No * 2; //number of total particles
    const double L = 11;//box length for internal materia

    double u[N], x[N], du[N];
    double rho = 100.0; //density and mass: constants
    double kappa = 1.0;//conductivity: constant
    double cv = 1.0;//heat compacity at consant volume
    double dx = L / ( Ni + 1 ), rc = 2.5 * dx;
    double m = dx * rho;
    double dt = 0.25 * dx * dx / kappa;

    //initial conditions
    for ( int i = 0; i < N; i ++ ) {
        x[i] = ( i - N / 2 ) * dx;
        u[i] = cv * ( i <= N / 2 ? 0 : 1 );
        if ( i >= N / 4 && i <= N * 3 / 4 ) u[i] += i % 2 * 0.5;
    }

    std::ofstream fout( "solution-sph.txt" );

    //time integration forward Euler scheme
    for ( int k = 0; k < 100; k++ ) {
        // push data to the other solver
        for ( int i : {3, 4, N - 5, N - 4} ) interface.push( "u", x[i], u[i] );
        interface.commit( k );

        sampler_gauss1d<double> gauss( rc, rc / 2 );
        chrono_sampler_exact1d  exact;
        for ( int i : {0, 1, N - 2, N - 1} ) u[i] = interface.fetch( "u", x[i], k, gauss, exact );

        //reset du
        for ( int i = 0; i < N; i ++ ) du[i] = 0.0;

        // N^2 brute-force pairwise evaluation
        for ( int i = 0; i < N - 1; i ++ ) {
            for ( int j = i + 1; j < N; j ++ ) {
                double r_ij  = x[j] - x[i];

                if ( r_ij <= rc ) {
                    double w_g = cubic_spline_gradient( r_ij, rc );
                    double du_ij = 2.0 * kappa * m / rho / rho * w_g * ( u[i] / cv - u[j] / cv );
                    du[i] += du_ij;
                    du[j] -= du_ij;
                }
            }
        }

        //only update the internal material points
        for ( int i = No; i < N - No; i ++ ) u[i] += du[i];

        if ( k % 10 == 0 ) {
            printf( "SPH step %d\n", k );
            for ( int i = No; i < N - No; i ++ ) fout << x[i] << '\t' << u[i] << std::endl;
            fout << std::endl;
        }
    }

    return 0;
}
