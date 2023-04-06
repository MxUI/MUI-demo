!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library Demo 10-1-5           *
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

  integer, parameter :: N = 110
  double precision :: u1(N), u2(N), k = 0.515, H = 1.0, rSearch = 1.0, tolerance = 8.0e-1
  type(c_ptr), target :: MUI_COMM_WORLD=c_null_ptr
  integer :: mui_ranks, mui_size, ierr, status, i, j, t, iter, pair_count
  type(c_ptr), target :: uniface_1d_f=c_null_ptr
  type(c_ptr), target :: mui_sampler_pseudo_nearest_neighbor_1d_f=c_null_ptr
  type(c_ptr), target :: mui_temporal_sampler_exact_1d_f=c_null_ptr
  type(c_ptr), target :: mui_algorithm_aitken_1d_f=c_null_ptr

  real(c_double), dimension (:), allocatable :: pp, value_init, u, v, tempValue
  character(len=1024) :: appname = "right"
  character(len=1024) :: uriheader = "mpi://"
  character(len=1024) :: uridomain = "/ifs"
  character(len=1024) :: name_fetch = "u"
  character(len=1024) :: name_push = "u0"
  character(len=1024) :: fileAddress = "results_right"
  character(len=1024) :: fileWriteName, fileWriteName2, fileWriteName3, fileWriteName4
  character(len=1024) :: uri1d, makedirMString
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
  print *, "{", trim(appname),"}: COMM_SIZE: ", mui_size, "COMM_RANK: ", mui_ranks

  do i=40,100,10
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
  do i = 40, 100, 10
      pp(i) = i
      value_init(i) = u1(i)
      pair_count = pair_count + 1
  end do

  ! MUI define spatial and temporal samplers
  call mui_create_sampler_pseudo_nearest_neighbor_1d_f(mui_sampler_pseudo_nearest_neighbor_1d_f, rSearch)
  call mui_create_temporal_sampler_exact_1d_f(mui_temporal_sampler_exact_1d_f, tolerance)
  call mui_create_algorithm_aitken_1d_f(mui_algorithm_fixed_relaxation_1d_f, DBLE(0.01), pp, value_init, pair_count)

  ! Create the file name
  write(fileWriteName, "(a, i0, a)") "results_right", mui_ranks, "/solution-right_AITKEN_0.csv"
  ! Open the file for writing
  open(unit=12, file=fileWriteName, action="write", status="replace")
  ! Write the header row
  write(12, "(a, a, a)") '"X"', ',','"u"'
  ! Write the data rows
  do i = 40, 100, 10
    write(12, "(f4.1, a, f4.2, a)")  i*H, ",", u(i), ","
  end do
  ! Close the file
  close(12)

  ! Create the file name
  write(fileWriteName3, "(a, i0, a)") "results_iteration_right", mui_ranks, "/solution-right_AITKEN_0.csv"
  ! Open the file for writing
  open(unit=14, file=fileWriteName3, action="write", status="replace")
  ! Write the header row
  write(14, "(a, a, a)") '"X"', ',','"u"'
  ! Write the data rows
  do i = 40, 100, 10
    write(14, "(f4.1, a, f4.2, a)")  i*H, ",", u(i), ","
  end do
  ! Close the file
  close(14)

  ! Begin time loops
  do t = 1, 10
      ! Begin iteration loops
      do iter = 1, 100

        write(*,*) "Right grid time ", t, " iteration ", iter

        ! MUI fetch points
        call mui_fetch_pseudo_nearest_neighbor_exact_fixed_relaxation_1d_pair_f(uniface_1d_f, &
                                       trim(name_fetch)//c_null_char, 40 * H, real(t, c_double), real(iter, c_double), &
                                        mui_sampler_pseudo_nearest_neighbor_1d_f,  &
                                        mui_temporal_sampler_exact_1d_f,mui_algorithm_fixed_relaxation_1d_f, u(40))

        if ((t .ge. 4) .and. (t .lt. 6)) then
            call mui_fetch_pseudo_nearest_neighbor_exact_fixed_relaxation_1d_pair_f(uniface_1d_f, &
                                           trim(name_fetch)//c_null_char, 42 * H, real(t, c_double), real(iter, c_double), &
                                            mui_sampler_pseudo_nearest_neighbor_1d_f,  &
                                            mui_temporal_sampler_exact_1d_f,mui_algorithm_fixed_relaxation_1d_f, u(42))
        end if

!        write(*,*) "Right under relaxation factor at t= ", t, " iter= ", iter, " is ", aitken_get_under_relaxation_factor(mui_algorithm_fixed_relaxation_1d_f,t,iter)
!        write(*,*) "Right residual L2 Norm at t= ", t, " iter= ", iter, " is ", aitken_get_residual_L2_Norm(mui_algorithm_fixed_relaxation_1d_f,t,iter)

        ! calculate 'interior' points
        do i = 50, 100, 10
            v(i) = u(i) + k / ( H * H ) * ( u(i - 10) + u(i + 10) - 2 * u(i) )
        end do

        ! calculate 'boundary' points
        v(N-10) = 0.0
        v(40) = u(40)

        if ((t .ge. 4) .and. (t .lt. 6)) then
            v(42) = u(42)
        end if

        ! MUI push points
        call mui_push_1d_f(uniface_1d_f, trim(name_push)//c_null_char, 60 * H, u(60))

        ! MUI commit
        call mui_commit_1d_pair_f(uniface_1d_f, real(t, c_double), real(iter, c_double))

        ! Swap u and v
        tempValue = u
        u = v
        v = tempValue

        ! Output
        write(fileWriteName2, "(a, i0, a, i0, a)") "results_iteration_right", mui_ranks,"/solution-right_AITKEN_", (((t-1)*100) + iter),".csv"
        open(unit=13, file=fileWriteName2, action="write", access="append")
        ! Write the header row
        write(13, "(a, a, a)") '"X"', ',','"u"'
        write(13, "(f4.1, a, f4.2, a)")  40*H, ",", u(40), ","
        if ((t .ge. 4) .and. (t .lt. 6)) then
            write(13, "(f4.1, a, f4.2, a)")  42*H, ",", u(42), ","
        end if
        do i = 50, 100, 10
            write(13, "(f4.1, a, f4.2, a)")  i*H, ",", u(i), ","
        end do
        close(unit=13)

      ! End iteration loop
      end do

      ! Output
      write(fileWriteName4, "(a, i0, a, i0, a)") "results_right", mui_ranks,"/solution-right_AITKEN_", t,".csv"
      open(unit=15, file=fileWriteName4, action="write", access="append")
      ! Write the header row
      write(15, "(a, a, a)") '"X"', ',','"u"'
      write(15, "(f4.1, a, f4.2, a)")  40*H, ",", u(40), ","
      if ((t .ge. 4) .and. (t .lt. 6)) then
        write(15, "(f4.1, a, f4.2, a)")  42*H, ",", u(42), ","
      end if
      do i = 50, 100, 10
        write(15, "(f4.1, a, f4.2, a)")  i*H, ",", u(i), ","
      end do
      close(unit=15)

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
