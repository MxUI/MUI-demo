// one dimension heat conductivity by SPH
// Author: Xin Bian
// Date  : 23:53:00, Sunday, June 19, 2016, from Munich.
//

#include <stdio.h>

double cubic_spline_gradient(double r, double rc);

main() {

  const int N = 19, Nint = 3;//number of inside/outside particles
  const int Ntotal= N + Nint * 2;//number of total particles
  const double L = 10;//box length for internal materia
  //rho du/dt = kappa Laplacian T

  double u[Ntotal], x[Ntotal], du[Ntotal];
  double rho = 1000.0; //density and mass: constants
  double kappa = 1.0;//conductivity: constant
  double cv = 1.0;//heat compacity at consant volume
  double dx = L / (N+1), rc = 2.5 * dx;
  double m = dx * rho;
  double T_left = 0.0, T_right = 1.0;
  double t_end = 0.5, dt = 0.25 * dx * dx / kappa;
  int step_end = 1000;

  //arbitrary initial conditions/boundary conditions
  //step_end = t_end/dt;
  for ( int i = 0; i< Ntotal; i ++ ) {
    x[i] = 0.0 + ( i + 1 - Nint ) * dx;
    if ( i <= Ntotal/2 ) {
      u[i] = cv * T_left;
    }
    else {
      u[i] = cv * T_right;
    }
  }

  //printf("step = %i\n", 0);
  for ( int i=0; i< Ntotal; i ++ ) {
    //printf("x[%i]=%f, u[%i]=%f\n", i, x[i], i, u[i]);
    printf("%f, %f\n", x[i], u[i]);
  }
  printf("\n");

  //time integration forward Euler scheme
  for ( int k = 0; k < step_end; k++ ) {

    //reset du
    for ( int i = 0; i < Ntotal; i ++ ) {
      du[i] = 0.0;
    }

    //straightforward pair-interactions
    //no cell or neighbor list for acceleration
    for ( int i = 0; i< Ntotal-1; i ++ ) {
      for ( int j = i+1; j < Ntotal; j ++ ) {
	//close enough
	double r_ij  = x[j] - x[i];

	if ( r_ij <= rc ) {

	  double w_g = cubic_spline_gradient(r_ij, rc);
	  double du_ij = 2.0*kappa*m /rho/rho * w_g * (u[i]/cv-u[j]/cv);

	  du[i] += du_ij;
	  du[j] -= du_ij;

	}// if ( r_ij <= rc )
      }//j
    }//i


    //only update the internal material points
    for ( int i = Nint; i < Ntotal-Nint; i ++ ) {
      u[i] += du[i];
    }

    //printf("step = %i\n", k+1);
    if ( k % 100 == 0 ) {
      for ( int i=0; i< Ntotal; i ++ ) {
	//printf("x[%i]=%f, u[%i]=%f\n", i, x[i], i, u[i]);
	printf("%f, %f\n", x[i], u[i]);
      }
      printf("\n");
    }
  }//k

}

double cubic_spline_gradient(double r, double rc) {
  const static double sigma  = 2.0 / 3.0; // for one diemension
  double h      = rc / 2.0;
  double coef_g = sigma / h / h;
  double s      = r / h;
  double w_g  = 0;
  if ( s < 1.0 ) w_g = -3.0*s + 9.0/4.0*s*s;
  else if ( s < 2.0 ) w_g = -3.0/4.0* (2-s) * (2-s);
  return w_g * coef_g;
}
