/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library Demo 6                *
*                                                                            *
* Copyright (C) 2019 S. Rolfo (*STFC Daresbury laboratory)                   *
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
 * @file DomainCoarse.cpp
 * @author Stefano Rolfo
 * @date May 16th 2019
 * @brief Simple programme to exchange an array using MUI with NN sampler
 * @brief The programme show how to start MUI in a funtion and COMMIT and FETCH in others
 * @brief Coarse Domain
 *
 * USAGE: sh run_case.sh
 */

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
  std::string interf_name="mpi://coarse/ifs";
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
  sampler_pseudo_nearest_neighbor2d<double> s1(1e-2);
  chrono_sampler_exact2d  s2;
  for(int i=0; i<nv; ++i){
    point2d loc( xx[i], yy[i] );
    ur[i] = interface->fetch("us",loc, it, s1, s2);
  }
}

/*
This file write the output in tecplot format forllowing a quadrilateral FE convention
*/
void output (int it, int nc, double *xx, double *yy,  double *us,  double *ur) 
{
  std::string root="solution_coarse";
  std::string format=".dat";
  std::string fname;
  // Convert to string and padding with 0
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









