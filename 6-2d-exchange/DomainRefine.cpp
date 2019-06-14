/*
 * DomainCoarse.cpp
 *
 *  Created on: May 16, 2019
 *      Author: S. Rolfo
 */

// Simple programme to exchange an array using MUI with NN sample
// The programme show how to start MUI in a funtion and COMMIT and FETCH in others
// Coarse Domain

#include "../mui/mui.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string> // for string class
#include <sstream>
#include <iomanip>

// include
typedef std::unique_ptr<mui::uniface2d> MuiUnifacePtr2d;
MuiUnifacePtr2d interface;

// Define functions
void init(int nv, double dx, double x0, double *xx, double *yy);
void compute_solution(int nv,int it,double tt,double *xx,double *yy,double *us);
void get_solution(int nv,int it,double tt,double *xx,double *yy,double *ur);
void output(int it,int nc, double *xx, double *yy, double *us, double *ur);

using namespace mui;


int main( int argc, char ** argv ) {

    // Definition of the paramters and declaration
    const static int    nc  = 9;
    const static int    nv  = nc+1;
    const static int    nv2 = nv*nv;
    const static int    nt  = 4;
    const static double dt  = 2.e-1;
    const static double lim[2]={-10.0, 10.0};
    double xx[nv2], yy[nv2], us[nv2], ur[nv2];  
    // Computing the delta space
    double dx = (lim[1]-lim[0]) / nc;
    // Init of the mesh 
    init (nv,dx,lim[0],xx,yy);
    //write of the results
    for(int it=1; it<=nt; ++it)
    {
      double tt = it*dt;
      compute_solution(nv2,it,tt,xx,yy,us);
      get_solution(nv2,it,tt,xx,yy,ur);
      if ((it % 2) == 0) output(it,nc,xx,yy,us,ur);
    }
    
    return 0;
}

void init (int nv, double dx, double x0, double *xx, double *yy) 
{
  // Init of MUI interface
  std::string interf_name="mpi://refine/ifs";
  interface.reset(new mui::uniface2d(interf_name));

  int i, j, id;

  for(j=0; j<nv; ++j) {
    for(i=0; i<nv; ++i){
      id = j*nv+i;
      xx[id] = x0+i*dx;
      yy[id] = x0+j*dx;
      //std::cout << xx[id] << " " <<yy[id] << std::endl;
    }
  }
}

void compute_solution(int nv,int it,double tt,double *xx,double *yy,double *us)
{
  for(int i=0; i<nv; ++i){
    us[i] = 1e1*tt*xx[i]*yy[i];
    point2d loc( xx[i], yy[i] );
    interface->push("us",loc,us[i]);
  }
  interface->commit( it );
}

void get_solution(int nv,int it,double tt,double *xx,double *yy,double *ur)
{
  sampler_nearest_neighbor2d<double> s1;
  chrono_sampler_exact2d  s2(1e-2);
  for(int i=0; i<nv; ++i){
    point2d loc( xx[i], yy[i] );
    ur[i] = interface->fetch("us",loc, it, s1, s2);
    //ur[i] = 1e1*tt*xx[i]*yy[i];
  }
}

/*
This file write the output in tecplot format forllowing a quadrilateral FE convention
*/
void output (int it, int nc, double *xx, double *yy,  double *us,  double *ur) 
{
  std::string root="solution_refine";
  std::string format=".dat";
  std::string fname;
  // Convert it in string padding with 0
  std::stringstream ss;
  ss << std::setw(3) << std::setfill('0') << it;
  std::string its = ss.str();
  // name of the file
  fname=root+its+format;
  std::ofstream fout(fname.c_str());
  // First line declaration variables
  fout << "Variables = x, y, uS, uR" << std::endl;
  fout << "  " << std::endl;
  fout << "Zone, N="<<(nc+1)*(nc+1)<<", E="<< nc*nc <<", F=FEBlock, ET=QUADRILATERAL" << std::endl;
  int nv  = (nc+1);
  int nv2 = nv*nv;
  for ( int i = 0; i < nv2; ++i ) fout << xx[i] << " ";
  fout << " \n";
  for ( int i = 0; i < nv2; ++i ) fout << yy[i] << " ";
  fout << " \n";
  for ( int i = 0; i < nv2; ++i ) fout << us[i] << " ";
  fout << " \n";
  for ( int i = 0; i < nv2; ++i ) fout << ur[i] << " ";
  fout << " \n";
  // write of the connectivity
  for ( int i = 1; i < nv2-nv; ++i ) 
  {
    if ((i % nv) != 0) fout <<i<<" "<<i+1<<" "<<i+1+nv<<" "<<i+nv<<" " ;
  }
  fout << " \n";
  fout.close();
}









