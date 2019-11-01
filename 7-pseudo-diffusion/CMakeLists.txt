# CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
 
PROJECT(pseudo_diffusion_RBF)

# <==> SET_INCLUDE_DIRs
FIND_PATH(MUIPATH mui.h)
IF(MUIPATH)
    INCLUDE_DIRECTORIES (${MUIPATH})
ELSE()
    MESSAGE( FATAL_ERROR "Please point the MUIPATH as 'export CMAKE_INCLUDE_PATH=/path/to/mui.h' to the include directory of your MUI code.")
    RETURN ()
ENDIF()

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

find_package (Eigen3 REQUIRED NO_MODULE)

if (TARGET Eigen3::Eigen)

    set(coarseSource 3D_pseudo_diffusion_coarse.cpp ${MUIPATH})
    set(fineSource 3D_pseudo_diffusion_fine.cpp ${MUIPATH})

    # <==> ADD_EXECUTABLE
    ADD_EXECUTABLE(3D_pseudo_diffusion_coarse ${coarseSource})
    target_compile_options (3D_pseudo_diffusion_coarse PRIVATE -std=c++11 -O3)
    target_compile_definitions (3D_pseudo_diffusion_coarse PRIVATE USE_RBF)
    target_link_libraries(3D_pseudo_diffusion_coarse ${MPI_LIBRARIES} Eigen3::Eigen)
    ADD_EXECUTABLE(3D_pseudo_diffusion_fine ${fineSource})
    target_compile_options (3D_pseudo_diffusion_fine PRIVATE -std=c++11 -O3)
    target_compile_definitions (3D_pseudo_diffusion_fine PRIVATE USE_RBF)
    target_link_libraries(3D_pseudo_diffusion_fine ${MPI_LIBRARIES} Eigen3::Eigen)
    ADD_EXECUTABLE(3D_pseudo_diffusion_standalone 3D_pseudo_diffusion_standalone.cpp)
    target_compile_options (3D_pseudo_diffusion_standalone PRIVATE -std=c++11 -O3)

    if(MPI_COMPILE_FLAGS)
      set_target_properties(3D_pseudo_diffusion_coarse PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
      set_target_properties(3D_pseudo_diffusion_fine PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif(MPI_COMPILE_FLAGS)

    if(MPI_LINK_FLAGS)
      set_target_properties(3D_pseudo_diffusion_coarse PROPERTIES
        LINK_FLAGS "${MPI_LINK_FLAGS}")
      set_target_properties(3D_pseudo_diffusion_fine PROPERTIES
        LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif(MPI_LINK_FLAGS)

ELSE()
    MESSAGE( FATAL_ERROR "Please point the path of Eigen3 as '-DEigen3_DIR=/path/to/your/eigen3/cmake/' to the include directory of your Eigen installation.")
    RETURN ()
ENDIF(TARGET Eigen3::Eigen)
