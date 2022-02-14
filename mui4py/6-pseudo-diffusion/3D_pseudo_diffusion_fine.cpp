/**
 *
 * @file 3D_pseudo_diffusion_fine.cpp
 * @author W.L.
 * @date 09 Oct 2019
 * @brief See README.md
 *
 */

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

/// Include MUI header file and configure file 
#include "mui.h"
#include "demo6_config.h"

int main(int argc, char **argv) {

  /// Create results folder
  mkdir("coupling_results", 0777);

  /// Create rbf matrix folder
  mkdir("rbfFineMatrix", 0777);

  /// Declare MPI common world with the scope of MUI
  MPI_Comm world = mui::mpi_split_by_app();

  /// Define the name of MUI domain
  std::string domain = "fineDomain";

  /// Define the name of MUI interfaces
  std::vector<std::string> interfaces;
  interfaces.emplace_back("interface2D01");
  interfaces.emplace_back("interface2D02");

  /// Declare MUI objects using MUI configure file
  auto ifs = mui::create_uniface < mui::demo6_config > (domain, interfaces);

  /// Declare MPI ranks and rank size
  int rank, size;
  MPI_Comm_rank(world, &rank);
  MPI_Comm_size(world, &size);

  /// Define the name of push/fetch values
  const char *name_push = "fineField";
  const char *name_fetch = "coarseField";

  /// Define the forget steps of MUI to reduce the memory
  int forgetSteps = 5;

  /// Define parameters of the RBF sampler
  double rSampler   = 0.4;                    // Define the search radius of the RBF sampler (radius size should be balanced to try and maintain)
  int basisFunc     = 0;                      // Specify RBF basis function 0-Gaussian; 1-WendlandC0; 2-WendlandC2; 3-WendlandC4; 4-WendlandC6
  bool conservative = false;                  // Enable conservative OR consistent RBF form
  bool polynomial   = true;                   // Enable/disable polynomial terms during RBF matrix creation
  bool smoothFunc   = false;                  // Enable/disable RBF smoothing function during matrix creation
  bool readMatrix   = false;                  // Enable/disable reading the matrix in from file
  bool writeMatrix  = true;                   // Enable/disable writing of the matrix (if not reading)
  std::string fileAddress("rbfFineMatrix");   // Output folder for the RBF matrix files
  double cutoff     = 1e-9;                   // Cut-off value for Gaussian RBF basis function
  double cgSolveTol = 1e-6;                   // Eigen Conjugate Gradient solver tolerance
  int cgMaxIter     = -1;                     // Eigen Conjugate Gradient solver maximum iterations (-1 = value determined by tolerance)
  int pouSize       = 50;                     // RBF Partition of Unity patch size

  /// Setup diffusion rate
  constexpr static double dr = 0.5;

  /// Setup time steps
  constexpr static int steps = 200;

  /// Setup output interval
  constexpr static int outputInterval = 1;

  /// Geometry info
  double x0 = 0.0; /// origin coordinate (x-axis direction) of the left part geometry
  double y0 = 0.0; /// origin coordinate (y-axis direction) of the left part geometry
  double z0 = 0.0; /// origin coordinate (z-axis direction) of the left part geometry
  double lx0 = 1.0; /// length (x-axis direction) of the left part geometry
  double ly0 = 1.0; /// length (y-axis direction) of the left part geometry
  double lz0 = 1.0; /// length (z-axis direction) of the left part geometry

  double x1 = 2.0; /// origin coordinate (x-axis direction) of the right part geometry
  double y1 = 0.0; /// origin coordinate (y-axis direction) of the right part geometry
  double z1 = 0.0; /// origin coordinate (z-axis direction) of the right part geometry
  double lx1 = 1.0; /// length (x-axis direction) of the right part geometry
  double ly1 = 1.0; /// length (y-axis direction) of the right part geometry
  double lz1 = 1.0; /// length (z-axis direction) of the right part geometry

  /// Domain discretization
  constexpr static int Nx = 11;     /// number of grid points in x axis per part
  constexpr static int Ny = 11;     /// number of grid points in y axis per part
  constexpr static int Nz = 11;     /// number of grid points in z axis per part
  constexpr static int Nt = 2 * Nx * Ny * Nz;  /// total number of points

  /// Declare points
  double points0[Nx][Ny][Nz][3], points1[Nx][Ny][Nz][3];

  /// Store point coordinates
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        points0[i][j][k][0] = x0 + (lx0 / (Nx - 1)) * i;
        points0[i][j][k][1] = y0 + (ly0 / (Ny - 1)) * j;
        points0[i][j][k][2] = z0 + (lz0 / (Nz - 1)) * k;

        points1[i][j][k][0] = x1 + (lx1 / (Nx - 1)) * i;
        points1[i][j][k][1] = y1 + (ly1 / (Ny - 1)) * j;
        points1[i][j][k][2] = z1 + (lz1 / (Nz - 1)) * k;
      }
    }
  }

  /// Generate initial pseudo scalar field
  double scalar_field0[Nx][Ny][Nz], scalar_field1[Nx][Ny][Nz];
  double tolerance = (lx0 / (Nx - 1)) * 0.5;
  double amplitude = 100.0;
  double x_center = x0 + lx0 / 2.0;
  double y_center = y0 + ly0 / 2.0;
  double z_center = z0 + lz0 / 2.0;
  double r_max = sqrt(
      pow((x0 - x_center), 2) + pow((y0 - y_center), 2)
          + pow((z0 - z_center), 2));

  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {

        if (points0[i][j][k][0] <= tolerance) {
          double r_ = sqrt(
              pow((points0[i][j][k][1] - y_center), 2)
                  + pow((points0[i][j][k][2] - z_center), 2));

          scalar_field0[i][j][k] = amplitude * (r_max - r_) / r_max;
        } else {
          scalar_field0[i][j][k] = 0.0;
        }
        scalar_field1[i][j][k] = 0.0;
      }
    }
  }

  /// Declare std::vector to store mui::point2d
  std::vector < mui::point2d > point2dvec;

  /// Store mui::point2d that located in the fetch interface
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        if ((points1[i][j][k][0] - x1) <= tolerance) {
          mui::point2d ptf(points1[i][j][k][1], points1[i][j][k][2]);
          point2dvec.push_back(ptf);
        }
      }
    }
  }

  /// Define and announce MUI send/receive span
  mui::geometry::box<mui::demo6_config> send_region( { y0, z0 },
      { (y0 + ly0), (z0 + lz0) });
  mui::geometry::box<mui::demo6_config> recv_region( { y1, z1 },
      { (y1 + ly1), (z1 + lz1) });
  ifs[0]->announce_send_span(0, steps, send_region);
  ifs[1]->announce_recv_span(0, steps, recv_region);

  /// Define spatial and temporal samplers
  mui::sampler_rbf<mui::demo6_config> spatial_sampler(rSampler, point2dvec,
        basisFunc, conservative, polynomial, smoothFunc, readMatrix, writeMatrix, fileAddress,
        cutoff, cgSolveTol, cgMaxIter, pouSize);
  mui::chrono_sampler_exact<mui::demo6_config> chrono_sampler;

  /// Commit ZERO step of MUI
  ifs[0]->commit(0);

  /// Output the initial pseudo scalar field
  std::ofstream outputFileLeft;
  std::string filenameL = "coupling_results/scalar_field_left_fine_0.csv";
  outputFileLeft.open(filenameL);
  outputFileLeft << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        outputFileLeft << points0[i][j][k][0] << "," << points0[i][j][k][1]
            << "," << points0[i][j][k][2] << "," << scalar_field0[i][j][k]
            << ", \n";
      }
    }
  }
  outputFileLeft.close();

  std::ofstream outputFileRight;
  std::string filenameR = "coupling_results/scalar_field_right_fine_0.csv";
  outputFileRight.open(filenameR);
  outputFileRight << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
  for (int k = 0; k < Nz; ++k) {
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        outputFileRight << points1[i][j][k][0] << "," << points1[i][j][k][1]
            << "," << points1[i][j][k][2] << "," << scalar_field1[i][j][k]
            << ", \n";
      }
    }
  }
  outputFileRight.close();

  /// Declare integrations of the pseudo scalar field at boundaries
  double intFaceLD1; // Integration of the pseudo scalar field at the left boundary of Domain 1
  double intFaceRD1; // Integration of the pseudo scalar field at the right boundary of Domain 1
  double intFaceLD3; // Integration of the pseudo scalar field at the left boundary of Domain 3
  double intFaceRD3; // Integration of the pseudo scalar field at the right boundary of Domain 3

  /// Define output files for boundary integrations
  std::ofstream outputIntegration;
  std::string outputIntegrationName =
      "coupling_results/faceIntegrationD1_D3.txt";
  outputIntegration.open(outputIntegrationName);
  outputIntegration
      << "\"t\",\"intFaceLD1\",\"intFaceRD1\",\"intFaceLD3\",\"intFaceRD3\"\n";
  outputIntegration.close();

  /// Begin time loops
  for (int t = 1; t <= steps; ++t) {

    printf("\n");
    printf("{Fine Domain} %d Step \n", t);

    /// Reset boundary integrations
    intFaceLD1 = 0.0;
    intFaceRD1 = 0.0;
    intFaceLD3 = 0.0;
    intFaceRD3 = 0.0;

    /// Loop over points of Domain 1
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

          /// Loop over 'internal' points of Domain 1
          if ((points0[i][j][k][0] - x0) > tolerance) {

            /// Calculate the diffusion of pseudo scalar field of Domain 1
            scalar_field0[i][j][k] += dr
                * (scalar_field0[(i - 1)][j][k] - scalar_field0[i][j][k]);

            /// Loop over right boundary points of Domain 1
            if (((x0 + lx0) - points0[i][j][k][0]) <= (tolerance)) {

              mui::point2d locp(points0[i][j][k][1], points0[i][j][k][2]);

              /// push data to the other solver
              ifs[0]->push(name_push, locp, scalar_field0[i][j][k]);

              /// Calculate the integration of right boundary points of Domain 1
              intFaceRD1 += scalar_field0[i][j][k];
            }

          } else { /// Loop over left boundary points of Domain 1

            /// Calculate the integration of left boundary points of Domain 1
            intFaceLD1 += scalar_field0[i][j][k];

          }
        }
      }
    }

    /// Commit 't' step of MUI
    int sent = ifs[0]->commit(t);

    /// Forget data from Zero step to 't - forgetSteps' step of MUI to save memory
    if ((t - forgetSteps) > 0) {
      ifs[0]->forget(t - forgetSteps);
    }

    /// Loop over points of Domain 3
    for (int k = 0; k < Nz; ++k) {
      for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {

          /// Loop over left boundary points of Domain 3
          if ((points1[i][j][k][0] - x1) <= tolerance) {

            mui::point2d locf(points1[i][j][k][1], points1[i][j][k][2]);

            /// Fetch data from the other solver
            scalar_field1[i][j][k] = ifs[1]->fetch(name_fetch, locf, t,
                spatial_sampler, chrono_sampler);

            /// Calculate the integration of left boundary points of Domain 3
            intFaceLD3 += scalar_field1[i][j][k];

          } else { /// Loop over 'internal' points of Domain 3

            /// Calculate the diffusion of pseudo scalar field of Domain 3
            scalar_field1[i][j][k] += dr
                * (scalar_field1[(i - 1)][j][k] - scalar_field1[i][j][k]);

            /// Loop over right boundary points of Domain 3
            if (((x1 + lx1) - points1[i][j][k][0]) <= (tolerance)) {

              /// Calculate the integration of right boundary points of Domain 3
              intFaceRD3 += scalar_field1[i][j][k];
            }
          }
        }
      }
    }

    /// Output the pseudo scalar field and the boundary integrations
    if ((t % outputInterval) == 0) {

      std::ofstream outputFileLeft;
      std::string filenameL = "coupling_results/scalar_field_left_fine_"
          + std::to_string(t) + ".csv";
      outputFileLeft.open(filenameL);
      outputFileLeft << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";

      for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
          for (int i = 0; i < Nx; ++i) {
            outputFileLeft << points0[i][j][k][0] << "," << points0[i][j][k][1]
                << "," << points0[i][j][k][2] << "," << scalar_field0[i][j][k]
                << ", \n";
          }
        }
      }
      outputFileLeft.close();

      std::ofstream outputFileRight;
      std::string filenameR = "coupling_results/scalar_field_right_fine_"
          + std::to_string(t) + ".csv";
      outputFileRight.open(filenameR);
      outputFileRight << "\"X\",\"Y\",\"Z\",\"scalar_field\"\n";
      for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
          for (int i = 0; i < Nx; ++i) {
            outputFileRight << points1[i][j][k][0] << "," << points1[i][j][k][1]
                << "," << points1[i][j][k][2] << "," << scalar_field1[i][j][k]
                << ", \n";
          }
        }
      }
      outputFileRight.close();

      outputIntegration.open(outputIntegrationName, std::ofstream::app);
      outputIntegration << std::to_string(t) << "," << intFaceLD1 << ","
          << intFaceRD1 << "," << intFaceLD3 << "," << intFaceRD3 << ", \n";
      outputIntegration.close();
    }
  }

  return 0;
}
