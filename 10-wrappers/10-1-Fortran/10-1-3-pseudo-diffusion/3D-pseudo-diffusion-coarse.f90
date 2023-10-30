!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library Demo 10-1-3           *
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
!Filename: 3D-pseudo-diffusion-coarse.f90
!Created: 05 April 2023
!Author:  W, Liu
!Description: Coarse (middle) domain of the 3D pseudo diffusion case

program main
    use iso_c_binding, only : c_ptr,c_null_char,c_double
    use iso_fortran_env, only : error_unit
    use mui_2d_f
    use mui_general_f

    implicit none

    ! MUI/MPI variables
    integer(c_int) :: num_interfaces = 2
    integer(c_int) :: outputInterval = 1
    integer(c_int) :: forgetSteps = 5
    integer(c_int) :: basisFunc = 1
    integer(c_int) :: conservative = 1
    integer(c_int) :: smoothFunc = 0
    integer(c_int) :: generateMatrix = 1
    integer(c_int) :: cgMaxIter = 500
    integer(c_int) :: preconditioner = 1
    integer(c_int) :: pouSize = 50
    integer(c_int) :: synchronised = 1
    integer(c_int) :: reset_log = 1
    real(c_double) :: rSampler   = 0.8_c_double
    real(c_double) :: cutoff   = 1.0e-9_c_double
    real(c_double) :: cgSolveTol = 1.0e-6_c_double
    real(c_double) :: tolerance_sampler=1e-37_c_double
    type(c_ptr), target :: spatial_sampler_rbf_2d=c_null_ptr
    type(c_ptr), target :: temporal_sampler_exact_2d=c_null_ptr
    integer :: mui_ranks, mui_size, ierr
    character(len=1024) :: numberSuffix
    character(:), allocatable, target :: interfaces(:)
    character(:), allocatable :: domain, name_fetchA, name_pushA
    type(c_ptr) :: MUI_COMM_WORLD

    ! Local parameters
    integer(c_int) :: Ntx = 9
    integer(c_int) :: Nty = 9
    integer(c_int) :: Ntz = 9
    integer(c_int) :: steps = 200
    integer :: i, j, k
    integer :: Nx, Ny, Nz, Nt, point_count, t
    real(c_double) :: dr = 0.5_c_double
    real(c_double), parameter :: lx = 1.0_c_double
    real(c_double), parameter :: ly = 1.0_c_double
    real(c_double), parameter :: x0 = 1.0_c_double
    real(c_double), parameter :: y0 = 0.0_c_double
    real(c_double) :: lzpr, lpz, lpz_half, z0pr, tolerance, intFaceLD2, intFaceRD2
    real(c_double) :: x, y, z
    real(c_double), dimension (:), allocatable :: point2dY, point2dZ
    real(c_double), dimension (:,:,:), allocatable :: scalar_field_A
    real(c_double), dimension (:,:,:,:), allocatable :: points

    character(len=32) :: makedirMString
    character(len=50) :: filename, filenameB, filenameC, filenameD
    character(len=:), allocatable :: fileAddress

    !Allociate memory based on number of interfaces
    ! Note: always allocate one more character for interfaces
    !   otherwise may cause "NULL interfaces" fatal error.
    allocate(character(len_trim("interface2D0")+2) :: interfaces(num_interfaces))
    !For multi-domain function, "uniface_pointers_2d" should be used to collect the array of
    ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
    allocate(uniface_pointers_2d(num_interfaces))

    !Call mui_mpi_split_by_app_f() function to init MPI
    call mui_mpi_split_by_app_f(MUI_COMM_WORLD)

    ! MUI/MPI get comm size & rank
    call mui_mpi_get_size_f(MUI_COMM_WORLD, mui_size)
    call mui_mpi_get_rank_f(MUI_COMM_WORLD, mui_ranks)

    if (mui_size > 2) then
        print *, "MPI Size larger than 2 does not supported yet."
        call MPI_Abort(MUI_COMM_WORLD, 1, ierr)
    end if

    Nx = Ntx                               ! number of grid points in x axis per MPI rank
    Ny = Nty                               ! number of grid points in y axis per MPI rank

    if (mui_ranks < mod(Ntz, mui_size)) then ! number of grid points in z axis per MPI rank
        Nz = Ntz/mui_size + 1
    else
        Nz = Ntz/mui_size
    end if

    Nt = Nx * Ny * Nz                      ! total number

    allocate (points(Nx,Ny,Nz,3), scalar_field_A(Nx,Ny,Nz))

    if(mui_size == 2) then
        if(mod(Ntz,2) == 0) then
            lpz = 1.0/dble(Ntz-1)
            lpz_half = lpz/2.0
            lzpr = (1.0/dble(mui_size)) - lpz_half ! length (z-axis direction) per MPI rank of the geometry per MPI rank
        else
            if(mui_ranks == 0) then
                lzpr = 1.0/dble(mui_size) ! length (z-axis direction) per MPI rank of the geometry per MPI rank
            else
                lpz = 1.0/dble(Ntz-1)
                lzpr = (1.0/dble(mui_size)) - lpz ! length (z-axis direction) per MPI rank of the geometry per MPI rank
            endif
        endif
    else if(mui_size == 1) then
        lzpr = 1.0
    endif

    if (mui_ranks == 1) then
        if (mod(Ntz, 2) == 0) then
            lpz = 1.0 / (dble(Ntz) - 1.0)
            lpz_half = lpz / 2.0
            z0pr = (1.0 / dble(mui_size)) + lpz_half ! origin coordinate (z-axis direction) per MPI rank of the geometry
        else
            lpz = 1.0 / (dble(Ntz) - 1.0)
            z0pr = (1.0 / dble(mui_size)) + lpz ! origin coordinate (z-axis direction) per MPI rank of the geometry
        end if
    else if (mui_ranks == 0) then
        z0pr = 0.0
    end if

    tolerance = (lx / (Nx - 1)) * 0.5

    ! Define rbf matrix folder
    write(makedirMString,'(a,i0)') "rbfCoarseMatrix", mui_ranks
    fileAddress = TRIM(makedirMString)

    ! Check if directory was created successfully
    if (ierr /= 0) then
        print *, "Error creating directory"
    else
        print *, "Directory created successfully"
    endif

    ! Define the name of MUI domain and fetch/push names
    domain = trim("coarseDomain")
    name_fetchA = trim("fineFieldA")
    name_pushA = trim("coarseFieldA")

    ! Create interface names
    do i = 1, num_interfaces
        !Generate character type of number suffix
        write (numberSuffix, "(I1)") i
        !Create and collect interface names
        interfaces(i) = trim("interface2D0") // trim(numberSuffix)
    end do

    !Create MUI interfaces. MUI interfaces will be collected by the "uniface_pointers_2d" after this subroutine
    call create_and_get_uniface_multi_2d_f(uniface_pointers_2d, domain, interfaces, num_interfaces)

    ! Store point coordinates
    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                points(i,j,k,1) = x0 + (lx / (Nx - 1)) * (i - 1)
                points(i,j,k,2) = y0 + (ly / (Ny - 1)) * (j - 1)
                points(i,j,k,3) = z0pr + (lzpr / (Nz - 1)) * (k - 1)
            end do
        end do
    end do

    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                scalar_field_A(i,j,k) = 0.0
            end do
        end do
    end do

    point_count = 0

    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                if ((points(i,j,k,1) - x0) <= tolerance) then
                    point_count = point_count + 1
                end if
            end do
        end do
    end do

    allocate (point2dY(point_count),point2dZ(point_count))

    point_count = 0

    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                if ((points(i,j,k,1) - x0) <= tolerance) then
                    point_count = point_count + 1
                    point2dY(point_count) = points(i,j,k,2)
                    point2dZ(point_count) = points(i,j,k,3)
                end if
            end do
        end do
    end do

    !Create spatial and temporal samplers for fetch operation
    call mui_create_sampler_rbf_2d_f(spatial_sampler_rbf_2d,rSampler,point2dY,point2dZ, &
                                        point_count,basisFunc,conservative,smoothFunc, &
                                        generateMatrix,TRIM(makedirMString),cutoff,cgSolveTol, &
                                        cgMaxIter,pouSize,preconditioner,MUI_COMM_WORLD)
    call mui_create_temporal_sampler_exact_2d_f(temporal_sampler_exact_2d, tolerance_sampler)

    call mui_commit_2d_f(uniface_pointers_2d(2)%ptr, DBLE(0))

    ! Output the initial pseudo scalar field
    write(filename, "(a, i0, a)") "coupling_results", mui_ranks, "/scalar_field_coarse_0.csv"
    ! Open the file for writing
    open(unit=11, file=filename, action="write", status="replace")
    write(11, "(a, a, a, a)") '"X","Y","Z","scalar_field_A"'
    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                x = points(i,j,k,1)
                y = points(i,j,k,2)
                z = points(i,j,k,3)
                write(11, "(f10.5, a, f10.5, a, f10.5, a, f10.5)") x, ',', y, ',', z, ',', scalar_field_A(i,j,k)
            end do
        end do
    end do
    ! Close the file
    close(11)

    ! Define output files for boundary integrations
    write(filenameB, "(a, i0, a)") "coupling_results", mui_ranks, "/faceIntegrationD2.txt"
    ! Open the file for writing
    open(unit=12, file=filenameB, action="write", status="replace")
    write(12, '(a)') '"t","intFaceLD2","intFaceRD2"'
    ! Close the file
    close(12)

    ! Begin time loops
    do t = 1, steps
        print *, ""
        print "(' {Coarse Domain} ', i0, ' Step')", t

        ! Reset boundary integrations
        intFaceLD2 = 0.0
        intFaceRD2 = 0.0

        ! Loop over points of Domain 2
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx

                    ! Loop over left boundary points of Domain 2
                    if ((points(i,j,k,1) - x0) <= tolerance) then

                        ! Fetch data from the other solvercall mui_push_2d_f(uniface_pointers_2d(i)%ptr, "position"//c_null_char, push_point_1, &
                        call mui_fetch_rbf_exact_2d_f(uniface_pointers_2d(1)%ptr,name_fetchA, &
                                                points(i,j,k,2),points(i,j,k,3),DBLE(t), &
                                                spatial_sampler_rbf_2d,temporal_sampler_exact_2d, &
                                                scalar_field_A(i,j,k))

                        ! Calculate the integration of left boundary points of Domain 2
                        intFaceLD2 = intFaceLD2 + scalar_field_A(i,j,k)

                    else ! Loop over 'internal' points of Domain 2

                        ! Calculate the diffusion of pseudo scalar field of Domain 2
                        scalar_field_A(i,j,k) = scalar_field_A(i,j,k) + dr * (scalar_field_A(i-1,j,k) - scalar_field_A(i,j,k))

                        ! Loop over right boundary points of Domain 2
                        if (((x0 + lx) - points(i,j,k,1)) <= tolerance) then

                            ! Push data to the other solver
                            call mui_push_2d_f(uniface_pointers_2d(2)%ptr, name_pushA, points(i,j,k,2), &
                                                points(i,j,k,3), scalar_field_A(i,j,k))

                            ! Calculate the integration of right boundary points of Domain 2
                            intFaceRD2 = intFaceRD2 + scalar_field_A(i,j,k)

                        end if
                    end if
                end do
            end do
        end do

        ! Commit 't' step of MUI
        call mui_commit_2d_f(uniface_pointers_2d(2)%ptr, DBLE(t))

        ! Output the pseudo scalar field and the boundary integrations
        if (mod(t, outputInterval) == 0) then
            write(filenameC, "(a, i0, a, i0, a)") "coupling_results", mui_ranks, "/scalar_field_coarse_", t, ".csv"
            ! Open the file for writing
            open(unit=13, file=filenameC, action="write", status="replace")
            write(13, "(A)") '"X","Y","Z","scalar_field_A"'
            do k = 1, Nz
                do j = 1, Ny
                    do i = 1, Nx
                        x = points(i,j,k,1)
                        y = points(i,j,k,2)
                        z = points(i,j,k,3)
                        write(13, "(f10.5, a, f10.5, a, f10.5, a, f10.5, a)") points(i,j,k,1), ',', &
                            points(i,j,k,2), ',', points(i,j,k,3), ',', scalar_field_A(i,j,k), ','
                    end do
                end do
            end do
            ! Close the file
            close(13)

            write(filenameD, "(a, i0, a)") "coupling_results", mui_ranks, "/faceIntegrationD2.txt"
            open(unit=14, file=filenameD, action="write", position='append')
            write(14, "(I0, a, f10.5, a, f10.5, a)") t, ',', intFaceLD2, ',', intFaceRD2, ','
            close(unit=14)

        end if

    ! End time loop
    end do

  !Destroy created 3D MUI objects
  call mui_destroy_sampler_rbf_2d_f(spatial_sampler_rbf_2d)
  call mui_destroy_temporal_sampler_exact_2d_f(temporal_sampler_exact_2d)
  ! Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  do i = 1, num_interfaces
     call mui_destroy_uniface_2d_f(uniface_pointers_2d(i)%ptr)
  end do
  !Release the memory on unifaces and interface names
  deallocate(uniface_pointers_2d)
  deallocate(interfaces)
  ! Deallocate arrays
  deallocate (points, scalar_field_A, point2dY, point2dZ, fileAddress)

end program main
