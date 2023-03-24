!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library                       *
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
!Description: Fortran demo to show smart send functions

program main
  use iso_c_binding, only : c_ptr,c_null_char,c_double
  use iso_fortran_env, only : error_unit
  use mui_3d_f
  use mui_general_f

  implicit none

  include 'mpif.h'

  !MUI/MPI variables
  real(c_double) :: zero=0.0_c_double
  integer :: reset_log = 1
  integer :: upper_forget = 5
  character(len=1024) :: URI
  type(c_ptr), target :: uniface_3t_f
  type(c_ptr), target :: spatial_sampler_pseudo_n2_linear_3t_f
  type(c_ptr), target :: temporal_sampler_exact_3t_f
  real(c_double) :: tolerance=8e-1_c_double
  integer :: memory_length = 5
  integer :: ierror, mui_ranks, mui_size
  type(c_ptr) :: MUI_COMM_WORLD

  ! Local parameters
  integer, parameter :: Nx        = 2 ! number of grid points in x axis
  integer, parameter :: Ny        = 2 ! number of grid points in y axis
  integer, parameter :: Nz        = 2 ! number of grid points in z axis
  character(len=1024) :: appname = "PUSHER_FETCHER_1"
  character(len=1024) :: uriheader = "mpi://"
  character(len=1024) :: uridomain = "/interface"
  character(len=1024) :: name_fetch = "displacement"
  character(len=1024) :: name_push = "pressure"
  real(c_double) :: rSearch    = 1.0 ! search radius
  integer :: Nt = Nx * Ny * Nz ! total time steps
  integer :: steps = 2 !  total time steps
  integer :: itersteps = 2 !  total iteration steps
  real(c_double) :: local_x0, local_y0, local_z0
  real(c_double) :: local_x1, local_y1, local_z1
  real(c_double) :: local_x2, local_y2, local_z2
  real(c_double) :: local_x3, local_y3, local_z3
  real(c_double), dimension (:,:,:,:), allocatable :: pp, pf
  real(c_double), dimension (:,:,:), allocatable :: pressure_push, pressure_fetch

  integer :: i, j, k, n, iter
  real(c_double) :: x, y, z

  allocate (pp(Nx,Ny,Nz,3), pf(Nx,Ny,Nz,3))  
  allocate (pressure_push(Nx,Ny,Nz), pressure_fetch(Nx,Ny,Nz))  

  !Call mui_mpi_split_by_app_f() function to init MPI
  call mui_mpi_split_by_app_f(MUI_COMM_WORLD)

  ! MUI set URL
  URI = trim(uriheader)//trim(appname)//trim(uridomain)
  print *, "{", trim(appname),"}: URI: ", trim(URI)

  ! MUI set uniface
  call mui_create_uniface_3t_f(uniface_3t_f, trim(URI)//c_null_char)

  ! MUI/MPI get comm size & rank
!  call MPI_COMM_SIZE(MUI_COMM_WORLD, mui_size, ierror)
!  call MPI_COMM_RANK(MUI_COMM_WORLD, mui_ranks, ierror)
!  print *, "{", trim(appname),"}: COMM_SIZE: ", mui_size, "COMM_RANK: ", mui_ranks

mui_size = 1
mui_ranks = 0

  ! define spatial and temporal samplers
  call mui_create_sampler_pseudo_n2_linear_3t_f(spatial_sampler_pseudo_n2_linear_3t_f, rSearch)
  call mui_create_temporal_sampler_exact_3t_f(temporal_sampler_exact_3t_f, tolerance)
!  call mui_set_forget_length_3t_f(uniface_3t_f(1),real(memory_length, c_double))
print *, "{", trim(appname),"}: HERE1 "
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
          pressure_push(i,j,k) = 32.2222
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
          pressure_fetch(i,j,k) = 11.11
       end do
    end do
  end do

  ! MUI annouce send/rcv span

  call mui_announce_send_span_3t_box_f(uniface_3t_f,local_x0,local_y0, &
    local_z0,local_x1,local_y1,local_z1,zero,real(steps, c_double),1)
  call mui_announce_recv_span_3t_box_f(uniface_3t_f,local_x2,local_y2, &
    local_z2,local_x3,local_y3,local_z3,zero,real(steps, c_double),1)

  ! Begin time loops
  do n = 1, steps

	do iter = 1, itersteps

		print *, "{", trim(appname),"}: ", n, " Step", iter, " Sub-iteration"

		! MUI push points
		do i = 1, Nx
		  do j = 1, Ny
			 do k = 1, Nz

				call mui_push_3t_f(uniface_3t_f, trim(name_push)//c_null_char, pp(i,j,k,1), pp(i,j,k,2), &
				pp(i,j,k,3), pressure_push(i,j,k))

			 end do
		  end do
		end do

		! MUI commit

		call mui_commit_3t_pair_f(uniface_3t_f, real(n, c_double), real(iter, c_double))

		! MUI fetch points
		do i = 1, Nx
		  do j = 1, Ny
			 do k = 1, Nz

				call mui_fetch_pseudo_n2_linear_exact_3t_pair_f(uniface_3t_f, &
												trim(name_fetch)//c_null_char, pf(i,j,k,1), pf(i,j,k,2), &
												pf(i,j,k,3), real(n, c_double), real(iter, c_double), &
												spatial_sampler_pseudo_n2_linear_3t_f, temporal_sampler_exact_3t_f, pressure_fetch(i,j,k))

			 end do
		  end do
		end do

		if ((n-upper_forget) .GT. 0) then
		  call mui_forget_upper_3t_pair_f(uniface_3t_f,real((n-upper_forget), c_double),real(itersteps, c_double),reset_log)
		end if

		! Print fetched values
		do i = 1, Nx
		  do j = 1, Ny
			 do k = 1, Nz
				print *, "{", trim(appname),"}: ", pressure_fetch(i,j,k)
			 end do
		  end do
		end do

	  ! End iteration loop
	  end do

  ! End time loop
  end do

  !Destroy created 3D MUI objects
  call mui_destroy_sampler_pseudo_nearest2_linear_3t_f(spatial_sampler_pseudo_n2_linear_3t_f)
  call mui_destroy_temporal_sampler_exact_3t_f(temporal_sampler_exact_3t_f)
  !Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  call mui_destroy_uniface_3t_f(uniface_3t_f)

  ! Deallocate arrays
  deallocate (pp, pf, pressure_push, pressure_fetch) 

end program main
