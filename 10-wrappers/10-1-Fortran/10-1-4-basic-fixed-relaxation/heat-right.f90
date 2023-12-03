!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library Demo 10-1-4           *
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
!Filename: heat-right.f90
!Created: 26 March 2023
!Author:  W, Liu
!Description: Right domain on basic Fixed Relaxation to demonstrate MUI coupling algorithm

program main
  use iso_c_binding, only : c_ptr,c_null_char,c_double,c_char, c_int32_t
  use iso_fortran_env, only : error_unit,output_unit, iostat_end
  use mui_1d_f
  use mui_general_f

!                            Left Domain                       Right Domain
! Coarse : +-------+-------+-------+-------o=======+=======o-------+-------+-------+-------+
!          0       1       2       3       4       5       6       7       8       9      10
! +: grid points
! o: interface points
! -: single domain zone
! =: overlapping zone
!

  implicit none

  include "mpif.h"

  integer, parameter :: N = 11
  double precision :: u1(N), u2(N), k = 0.515, H = 1.0, rSearch = 1.0, tolerance = 8.0e-1
  integer(c_int) ::  MUI_COMM_WORLD
  integer :: mui_ranks, mui_size, ierr, status, i, j, iter, pair_count, my_rank, num_procs, total_rank, total_procs
  type(c_ptr), target :: uniface_1d_f=c_null_ptr
  type(c_ptr), target :: mui_sampler_pseudo_nearest_neighbor_1d_f=c_null_ptr
  type(c_ptr), target :: mui_temporal_sampler_exact_1d_f=c_null_ptr
  type(c_ptr), target :: mui_algorithm_fixed_relaxation_1d_f=c_null_ptr

  real(c_double), dimension (:), allocatable :: pp, value_init, u, v, tempValue
  real(c_double) :: underRelax,resL2Norm
  character(len=1024) :: appname = "right"
  character(len=1024) :: uriheader = "mpi://"
  character(len=1024) :: uridomain = "/ifs"
  character(len=1024) :: name_fetch = "u"
  character(len=1024) :: name_push = "u0"
  character(len=1024) :: fileAddress = "results_right"
  character(len=1024) :: fileWriteName, fileWriteName2, uri1d, makedirMString
  logical :: fileExists

  allocate (pp(N), value_init(N), u(N), v(N), tempValue(N))

  ! MUI set URL
  uri1d = trim(uriheader)//trim(appname)//trim(uridomain)
  print *, "{", trim(appname),"}: URI: ", trim(uri1d)

  ! MUI set uniface
  call mui_create_uniface_1d_f(uniface_1d_f, trim(uri1d)//c_null_char)

  !Call mui_mpi_split_by_app_f() function to init MPI
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

  do i = 4, 10
    u1(i) = 0.0
  end do

  ! MUI/MPI get comm size & rank
  call mui_mpi_get_size_f(MUI_COMM_WORLD, mui_size)
  call mui_mpi_get_rank_f(MUI_COMM_WORLD, mui_ranks)
  print *, "{", trim(appname),"}: COMM_SIZE: ", mui_size, "COMM_RANK: ", mui_ranks

  ! Convert integer to string
  write(makedirMString,'(a,i0)') "results_right", mui_ranks
  ! Create directory
  call system('mkdir -m 777 '//trim(makedirMString), ierr)

  ! Check if directory was created successfully
  if (ierr /= 0) then
    print *, "Error creating directory"
  else
    print *, "Directory created successfully"
  endif

  ! Set push points & value
  pair_count = 0
  do i = 4, 10
      pp(i) = i
      value_init(i) = u1(i)
      pair_count = pair_count + 1
  end do

  ! MUI define spatial and temporal samplers
  call mui_create_sampler_pseudo_nearest_neighbor_1d_f(mui_sampler_pseudo_nearest_neighbor_1d_f, rSearch)
  call mui_create_temporal_sampler_exact_1d_f(mui_temporal_sampler_exact_1d_f, tolerance)
  call mui_create_algorithm_fixed_relaxation_1d_f(mui_algorithm_fixed_relaxation_1d_f, DBLE(0.01), & 
                                                  MUI_COMM_WORLD, pp, value_init, pair_count)

  ! Create the file name
  write(fileWriteName, "(a, i0, a)") "results_right", mui_ranks, "/solution-right_FR_0.csv"
  ! Open the file for writing
  open(unit=12, file=fileWriteName, action="write", status="replace")
  ! Write the header row
  write(12, "(a, a, a)") '"X"', ',','"u"'
  ! Write the data rows
  do i = 4, 10
    write(12, "(f4.1, a, f4.2, a)")  i*H, ",", u(i), ","
  end do
  ! Close the file
  close(12)

  ! Begin time loops
  do iter = 1, 1000

    write(*,*) "Right grid iteration ", iter

    ! MUI fetch points
    call mui_fetch_pseudo_n2_exact_fixed_relaxation_1d_f(uniface_1d_f, &
                                   trim(name_fetch)//c_null_char, 4 * H, real(iter, c_double), &
                                    mui_sampler_pseudo_nearest_neighbor_1d_f,  &
                                    mui_temporal_sampler_exact_1d_f,mui_algorithm_fixed_relaxation_1d_f, u(4))

    call mui_fixed_relaxation_get_under_relaxation_factor_1d_f(mui_algorithm_fixed_relaxation_1d_f,  &
                                        real(iter, c_double), underRelax)
    call mui_fixed_relaxation_get_residual_1d_f(mui_algorithm_fixed_relaxation_1d_f,  &
                                        real(iter, c_double), resL2Norm)

    write(*,*) "Right under relaxation factor at iter= ", iter, " is ", underRelax
    write(*,*) "Right residual L2 Norm at iter= ", iter, " is ", resL2Norm

    ! calculate 'interior' points
    do i = 5, 10
        v(i) = u(i) + k / ( H * H ) * ( u(i - 1) + u(i + 1) - 2 * u(i) )
    end do

    ! calculate 'boundary' points
    v(N-1) = 0.0
    v(4) = u(4)

    ! MUI push points
    call mui_push_1d_f(uniface_1d_f, trim(name_push)//c_null_char, 6 * H, u(6))

    ! MUI commit
    call mui_commit_1d_f(uniface_1d_f, real(iter, c_double))

    ! Swap u and v
    tempValue = u
    u = v
    v = tempValue

    ! Output
    write(fileWriteName2, "(a, i0, a, i0, a)") "results_right", mui_ranks,"/solution-right_FR_0", iter,".csv"
    open(unit=13, file=fileWriteName2, action="write", access="append")
    ! Write the header row
    write(13, "(a, a, a)") '"X"', ',','"u"'
    do i = 4, 10
        write(13, "(f4.1, a, f4.2, a)")  i*H, ",", u(i), ","
    end do
    close(unit=13)

  ! End time loop
  end do

  !Destroy created 3D MUI objects
  call mui_destroy_sampler_pseudo_nearest_neighbor_1d_f(mui_sampler_pseudo_nearest_neighbor_1d_f)
  call mui_destroy_temporal_sampler_exact_1d_f(mui_temporal_sampler_exact_1d_f)
  !Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  call mui_destroy_uniface_1d_f(uniface_1d_f)
  ! Deallocate arrays
  deallocate (pp, value_init)

end program main
