!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library Demo 10-1-2           *
!*                                                                            *
!* Copyright (C) 2023 W. Liu                                                  *
!*                                                                            *
!* This software is jointly licensed under the Apache License, Version 2.0    *
!* and the GNU General Public License version 3, you may use it according     *
!* to either.                                                                 *
!*                                                                            *
!* ** Apache License, version 2.0 **                                          *
!*                                                                            *
!* Licensed under the Apache License, Version 2.0 (the "License");            *
!* you may not use this file except in compliance with the License.           *
!* You may obtain a copy of the License at                                    *
!*                                                                            *
!* http://www.apache.org/licenses/LICENSE-2.0                                 *
!*                                                                            *
!* Unless required by applicable law or agreed to in writing, software        *
!* distributed under the License is distributed on an "AS IS" BASIS,          *
!* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
!* See the License for the specific language governing permissions and        *
!* limitations under the License.                                             *
!*                                                                            *
!* ** GNU General Public License, version 3 **                                *
!*                                                                            *
!* This program is free software: you can redistribute it and/or modify       *
!* it under the terms of the GNU General Public License as published by       *
!* the Free Software Foundation, either version 3 of the License, or          *
!* (at your option) any later version.                                        *
!*                                                                            *
!* This program is distributed in the hope that it will be useful,            *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
!* GNU General Public License for more details.                               *
!*                                                                            *
!* You should have received a copy of the GNU General Public License          *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
!******************************************************************************
!
!** File Details **
!
!Filename: code1.f90
!Created: 23 March 2023
!Author:  W, Liu
!Description: Fortran demo to show smart send functions with dual time types

program main
  use iso_c_binding, only : c_ptr,c_null_char,c_double
  use iso_fortran_env, only : error_unit
  use mui_3d_f
  use mui_general_f

  implicit none

  include "mpif.h"

  ! MUI/MPI variables
  real(c_double) :: zero=0.0_c_double
  real(c_double) :: tolerance=8e-1_c_double
  integer :: reset_log = 1
  integer :: upper_forget = 5
  integer :: synchronised = 1
  integer :: mui_ranks, mui_size, ierr, my_rank, num_procs, total_rank, total_procs
  character(len=1024) :: URI
  integer(c_int) :: MUI_COMM_WORLD
  type(c_ptr), target :: uniface_3d_f=c_null_ptr
  type(c_ptr), target :: spatial_sampler_pseudo_n2_linear_3d_f=c_null_ptr
  type(c_ptr), target :: temporal_sampler_exact_3d_f=c_null_ptr

  ! Local parameters
  integer, parameter :: Nx = 2 ! number of grid points in x axis
  integer, parameter :: Ny = 2 ! number of grid points in y axis
  integer, parameter :: Nz = 2 ! number of grid points in z axis
  integer, parameter :: steps = 2 !  total time steps
  integer, parameter :: itersteps = 2 !  total iteration steps
  real(c_double), parameter :: rSearch = 1.0 ! search radius
  character(len=1024) :: appname = "PUSHER_FETCHER_1"
  character(len=1024) :: uriheader = "mpi://"
  character(len=1024) :: uridomain = "/interface"
  character(len=1024) :: name_fetch = "code2Value"
  character(len=1024) :: name_push = "code1Value"

  integer :: Nt = Nx * Ny * Nz
  integer :: i, j, k, n, iter
  real(c_double) :: x, y, z
  real(c_double) :: local_x0, local_y0, local_z0
  real(c_double) :: local_x1, local_y1, local_z1
  real(c_double) :: local_x2, local_y2, local_z2
  real(c_double) :: local_x3, local_y3, local_z3
  real(c_double), dimension (:,:,:,:), allocatable :: pp, pf
  real(c_double), dimension (:,:,:), allocatable :: value_push, value_fetch

  allocate (pp(Nx,Ny,Nz,3), pf(Nx,Ny,Nz,3))  
  allocate (value_push(Nx,Ny,Nz), value_fetch(Nx,Ny,Nz))

  ! Call mui_mpi_split_by_app_f() function to init MPI
  call mui_mpi_split_by_app_f(MUI_COMM_WORLD)

  ! MUI/MPI get comm size & rank
  call mui_mpi_get_size_f(MUI_COMM_WORLD, mui_size)
  call mui_mpi_get_rank_f(MUI_COMM_WORLD, mui_ranks)

  call MPI_COMM_RANK(MUI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MUI_COMM_WORLD, num_procs, ierr)

  call MPI_COMM_RANK(MPI_COMM_WORLD, total_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, total_procs, ierr)

  print *, "{", trim(appname),"}: MPI Size from MUI func: ", mui_size, &
           " MPI Size from MPI func with MUI_COMM_WORLD: ", num_procs, &
           " MPI Size from MPI func with MPI_COMM_WORLD: ", total_procs
  print *, "{", trim(appname),"}: MPI Rank from MUI func: ", mui_ranks, &
           " MPI Rank from MPI func with MUI_COMM_WORLD: ", my_rank, &
           " MPI Rank from MPI func with MPI_COMM_WORLD: ", total_rank

  ! MUI set URL
  URI = trim(uriheader)//trim(appname)//trim(uridomain)
  print *, "{", trim(appname),"}: URI: ", trim(URI)

  ! MUI set uniface
  call mui_create_uniface_3d_f(uniface_3d_f, trim(URI)//c_null_char)

  ! Define spatial and temporal samplers
  call mui_create_sampler_pseudo_n2_linear_3d_f(spatial_sampler_pseudo_n2_linear_3d_f, rSearch)
  call mui_create_temporal_sampler_exact_3d_f(temporal_sampler_exact_3d_f, tolerance)

  ! Define bounding box
  local_x0 = 6.0 !  local push origin box start
  local_y0 = 0.0
  local_z0 = 0.0 + ((2.0/real(mui_size, c_double))*real(mui_ranks, c_double))
  local_x1 = 7.0 !  local push origin box end
  local_y1 = 1.0
  local_z1 = 0.0 + ((2.0/real(mui_size, c_double))*(real(mui_ranks, c_double)+1))
  local_x2 = 4.0 !  local fetch origin box start
  local_y2 = 0.0
  local_z2 = 0.0 + ((2.0/real(mui_size, c_double))*real(mui_ranks, c_double))
  local_x3 = 5.0 !  local fetch origin box end
  local_y3 = 1.0
  local_z3 = 0.0 + ((2.0/real(mui_size, c_double))*(real(mui_ranks, c_double)+1))

  ! Set push points & value
  do i = 1, Nx
    do j = 1, Ny
       do k = 1, Nz
          x = local_x0 + (((local_x1 - local_x0)/real(Nx, c_double))*real((i-1), c_double))
          y = local_y0 + (((local_y1 - local_y0)/real(Ny, c_double))*real((j-1), c_double))
          z = local_z0 + (((local_z1 - local_z0)/real(Nz, c_double))*real((k-1), c_double))
          pp(i,j,k,1) = x
          pp(i,j,k,2) = y
          pp(i,j,k,3) = z
          value_push(i,j,k) = 32.2222
       end do
    end do
  end do

  ! Initialise fetch points & value
  do i = 1, Nx
    do j = 1, Ny
       do k = 1, Nz
          x = local_x2 + (((local_x3 - local_x2)/real(Nx, c_double))*real((i-1), c_double))
          y = local_y2 + (((local_y3 - local_y2)/real(Ny, c_double))*real((j-1), c_double))
          z = local_z2 + (((local_z3 - local_z2)/real(Nz, c_double))*real((k-1), c_double))
          pf(i,j,k,1) = x
          pf(i,j,k,2) = y
          pf(i,j,k,3) = z
          value_fetch(i,j,k) = 11.11
       end do
    end do
  end do

  ! MUI annouce send/rcv span
  call mui_announce_send_span_3d_box_f(uniface_3d_f,local_x0,local_y0, &
    local_z0,local_x1,local_y1,local_z1,zero,real(steps, c_double),synchronised)
  call mui_announce_recv_span_3d_box_f(uniface_3d_f,local_x2,local_y2, &
    local_z2,local_x3,local_y3,local_z3,zero,real(steps, c_double),synchronised)

  ! Begin time loops
  do n = 1, steps
	do iter = 1, itersteps

		print *, "{", trim(appname),"}: ", n, " Step", iter, " Sub-iteration"

		! MUI push points
		do i = 1, Nx
		  do j = 1, Ny
			 do k = 1, Nz
				call mui_push_3d_f(uniface_3d_f, trim(name_push)//c_null_char, pp(i,j,k,1), pp(i,j,k,2), &
				pp(i,j,k,3), value_push(i,j,k))
			 end do
		  end do
		end do

		! MUI commit
		call mui_commit_3d_pair_f(uniface_3d_f, real(n, c_double), real(iter, c_double))

		! MUI fetch points
		do i = 1, Nx
		  do j = 1, Ny
			 do k = 1, Nz
				call mui_fetch_pseudo_n2_linear_exact_3d_pair_f(uniface_3d_f, &
												trim(name_fetch)//c_null_char, pf(i,j,k,1), pf(i,j,k,2), &
												pf(i,j,k,3), real(n, c_double), real(iter, c_double), &
												spatial_sampler_pseudo_n2_linear_3d_f, temporal_sampler_exact_3d_f, value_fetch(i,j,k))
			 end do
		  end do
		end do

		if ((n-upper_forget) .GT. 0) then
		  call mui_forget_upper_3d_pair_f(uniface_3d_f,real((n-upper_forget), c_double),real(itersteps, c_double),reset_log)
		end if

		! Print fetched values
		do i = 1, Nx
		  do j = 1, Ny
			 do k = 1, Nz
				print *, "{", trim(appname),"}: ", value_fetch(i,j,k)
			 end do
		  end do
		end do

	  ! End iteration loop
	  end do
  ! End time loop
  end do

  !Destroy created 3D MUI objects
  call mui_destroy_sampler_pseudo_n2_linear_3d_f(spatial_sampler_pseudo_n2_linear_3d_f)
  call mui_destroy_temporal_sampler_exact_3d_f(temporal_sampler_exact_3d_f)
  !Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  call mui_destroy_uniface_3d_f(uniface_3d_f)
  ! Deallocate arrays
  deallocate (pp, pf, value_push, value_fetch)

end program main
