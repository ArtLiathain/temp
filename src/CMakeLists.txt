#
#  This file is part of the FFEA simulation package
#
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file.
#
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
#
#  To help us fund FFEA development, we humbly ask that you cite
#  the research papers on the package.
#

include_directories("${PROJECT_SOURCE_DIR}/include")

# Tracking headers as source files, ensures the target is seen to depend on them and will rebuild if they change
set(FFEA_LIB_INCLUDE
    ${PROJECT_SOURCE_DIR}/include/FFEA_return_codes.h
    ${PROJECT_SOURCE_DIR}/include/FFEA_input_reader.h
    ${PROJECT_SOURCE_DIR}/include/BEM_Poisson_Boltzmann.h
    ${PROJECT_SOURCE_DIR}/include/mesh_node.h
    ${PROJECT_SOURCE_DIR}/include/SimulationParams.h
    ${PROJECT_SOURCE_DIR}/include/tetra_element_linear.h
    ${PROJECT_SOURCE_DIR}/include/SparseSubstitutionSolver.h
    ${PROJECT_SOURCE_DIR}/include/Face.h
    ${PROJECT_SOURCE_DIR}/include/BiCGSTAB_solver.h
    ${PROJECT_SOURCE_DIR}/include/GaussianQuadrature_1d.h
    ${PROJECT_SOURCE_DIR}/include/GaussianQuadrature_tri.h
    ${PROJECT_SOURCE_DIR}/include/MassLumpedSolver.h
    ${PROJECT_SOURCE_DIR}/include/MassMatrixLinear.h
    ${PROJECT_SOURCE_DIR}/include/MassMatrixQuadratic.h
    ${PROJECT_SOURCE_DIR}/include/MatrixFixedSparsityPattern.h
    ${PROJECT_SOURCE_DIR}/include/LJ_matrix.h
    ${PROJECT_SOURCE_DIR}/include/BindingSite.h
    ${PROJECT_SOURCE_DIR}/include/mat_vec_fns_II.h
    ${PROJECT_SOURCE_DIR}/include/NearestNeighbourLinkedListCube.h
    ${PROJECT_SOURCE_DIR}/include/PoissonMatrixQuadratic.h
    ${PROJECT_SOURCE_DIR}/include/SparseMatrixTypes.h
    ${PROJECT_SOURCE_DIR}/include/SparseMatrixUnknownPattern.h
    ${PROJECT_SOURCE_DIR}/include/SparsityPattern.h
    ${PROJECT_SOURCE_DIR}/include/Spring.h
    ${PROJECT_SOURCE_DIR}/include/VdW_solver.h
    ${PROJECT_SOURCE_DIR}/include/Steric_solver.h
    ${PROJECT_SOURCE_DIR}/include/LJSteric_solver.h
    ${PROJECT_SOURCE_DIR}/include/GenSoftSSINT_solver.h
    ${PROJECT_SOURCE_DIR}/include/KineticState.h
    ${PROJECT_SOURCE_DIR}/include/PreComp_solver.h
    ${PROJECT_SOURCE_DIR}/include/CheckTetrahedraOverlap.h
    ${PROJECT_SOURCE_DIR}/include/VolumeIntersection.h
    ${PROJECT_SOURCE_DIR}/include/RngStream.h
    ${PROJECT_SOURCE_DIR}/include/rod_blob_interface.h
    ${PROJECT_SOURCE_DIR}/include/ffea_test.h
)
set(FFEA_INCLUDE
    ${PROJECT_SOURCE_DIR}/include/CG_solver.h
    ${PROJECT_SOURCE_DIR}/include/ConjugateGradientSolver.h
    ${PROJECT_SOURCE_DIR}/include/FFEA_user_info.h
    ${PROJECT_SOURCE_DIR}/include/mat_vec_fns.h
    ${PROJECT_SOURCE_DIR}/include/NoMassCGSolver.h
    ${PROJECT_SOURCE_DIR}/include/SparseMatrixFixedPattern.h
    ${PROJECT_SOURCE_DIR}/include/World.h
    ${PROJECT_SOURCE_DIR}/include/Blob.h
    ${PROJECT_SOURCE_DIR}/include/rod_math_v9.h
    ${PROJECT_SOURCE_DIR}/include/rod_interactions.h
    ${PROJECT_SOURCE_DIR}/include/rod_structure.h
)
# These headers don't have a matching source file, so currently assumed part of ffea_lib
set(UNKNOWN_INCLUDE
    ${PROJECT_SOURCE_DIR}/include/BlobLite.h
    ${PROJECT_SOURCE_DIR}/include/ConnectivityTypes.h
    ${PROJECT_SOURCE_DIR}/include/LinkedListCube.h
    ${PROJECT_SOURCE_DIR}/include/NodeMap.h
    ${PROJECT_SOURCE_DIR}/include/SecondOrderFunctions.h
    ${PROJECT_SOURCE_DIR}/include/Solver.h
    ${PROJECT_SOURCE_DIR}/include/dimensions.h
    ${PROJECT_SOURCE_DIR}/include/mat_vec_types.h
    ${PROJECT_BINARY_DIR}/include/FFEA_version.h # Dynamically created from template at CMake configure
)

set(FFEA_LIB_SRC
    ${PROJECT_SOURCE_DIR}/src/FFEA_return_codes.cpp
    ${PROJECT_SOURCE_DIR}/src/FFEA_input_reader.cpp
    ${PROJECT_SOURCE_DIR}/src/BEM_Poisson_Boltzmann.cpp
    ${PROJECT_SOURCE_DIR}/src/mesh_node.cpp
    ${PROJECT_SOURCE_DIR}/src/SimulationParams.cpp
    ${PROJECT_SOURCE_DIR}/src/tetra_element_linear.cpp
    ${PROJECT_SOURCE_DIR}/src/SparseSubstitutionSolver.cpp
    ${PROJECT_SOURCE_DIR}/src/Face.cpp
    ${PROJECT_SOURCE_DIR}/src/BiCGSTAB_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/GaussianQuadrature_1d.cpp
    ${PROJECT_SOURCE_DIR}/src/GaussianQuadrature_tri.cpp
    ${PROJECT_SOURCE_DIR}/src/MassLumpedSolver.cpp
    ${PROJECT_SOURCE_DIR}/src/MassMatrixLinear.cpp
    ${PROJECT_SOURCE_DIR}/src/MassMatrixQuadratic.cpp
    ${PROJECT_SOURCE_DIR}/src/MatrixFixedSparsityPattern.cpp
    ${PROJECT_SOURCE_DIR}/src/LJ_matrix.cpp
    ${PROJECT_SOURCE_DIR}/src/BindingSite.cpp
    ${PROJECT_SOURCE_DIR}/src/mat_vec_fns_II.cpp
    ${PROJECT_SOURCE_DIR}/src/NearestNeighbourLinkedListCube.cpp
    ${PROJECT_SOURCE_DIR}/src/PoissonMatrixQuadratic.cpp
    ${PROJECT_SOURCE_DIR}/src/SparseMatrixTypes.cpp
    ${PROJECT_SOURCE_DIR}/src/SparseMatrixUnknownPattern.cpp
    ${PROJECT_SOURCE_DIR}/src/SparsityPattern.cpp
    ${PROJECT_SOURCE_DIR}/src/Spring.cpp
    ${PROJECT_SOURCE_DIR}/src/VdW_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/Steric_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/LJSteric_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/GenSoftSSINT_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/KineticState.cpp
    ${PROJECT_SOURCE_DIR}/src/PreComp_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/CheckTetrahedraOverlap.cpp
    ${PROJECT_SOURCE_DIR}/src/VolumeIntersection.cpp
    ${PROJECT_SOURCE_DIR}/src/RngStream.cpp
    #${PROJECT_SOURCE_DIR}/src/rod_math_v9.cpp
    #${PROJECT_SOURCE_DIR}/src/rod_interactions.cpp
    #${PROJECT_SOURCE_DIR}/src/rod_structure.cpp
    ${PROJECT_SOURCE_DIR}/src/rod_blob_interface.cpp
    ${PROJECT_SOURCE_DIR}/src/ffea_test.cpp
)
set(FFEA_SRC
    ${PROJECT_SOURCE_DIR}/src/ffea.cpp
    ${PROJECT_SOURCE_DIR}/src/CG_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/ConjugateGradientSolver.cpp
    ${PROJECT_SOURCE_DIR}/src/FFEA_user_info.cpp
    ${PROJECT_SOURCE_DIR}/src/mat_vec_fns.cpp
    ${PROJECT_SOURCE_DIR}/src/NoMassCGSolver.cpp
    ${PROJECT_SOURCE_DIR}/src/SparseMatrixFixedPattern.cpp
    ${PROJECT_SOURCE_DIR}/src/World.cpp
    ${PROJECT_SOURCE_DIR}/src/Blob.cpp
    ${PROJECT_SOURCE_DIR}/src/rod_math_v9.cpp
    ${PROJECT_SOURCE_DIR}/src/rod_interactions.cpp
    ${PROJECT_SOURCE_DIR}/src/rod_structure.cpp
)

add_library(ffea_lib ${FFEA_LIB_SRC} ${FFEA_LIB_INCLUDE} ${UNKNOWN_INCLUDE})
target_compile_features(ffea_lib PUBLIC cxx_std_17)
# target_link_libraries(ffea_lib PRIVATE stdc++fs) # Core to stdlib as of C++17
# Make the include dir available to dependencies
set(FFEA_INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include
)
set_target_properties(ffea_lib PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFEA_INCLUDE_DIRS}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ffea ${FFEA_SRC} ${FFEA_INCLUDE})
target_link_libraries(ffea PRIVATE ffea_lib)


if(USE_OPENMP)
    target_compile_definitions(ffea PRIVATE FFEA_PARALLEL_WITHIN_BLOB)
    add_executable(ffea_mb
        ffea.cpp
        CG_solver.cpp
        ConjugateGradientSolver.cpp
        FFEA_user_info.cpp
        mat_vec_fns.cpp
        NoMassCGSolver.cpp
        SparseMatrixFixedPattern.cpp
        World.cpp
        Blob.cpp
        rod_math_v9.cpp
        rod_interactions.cpp
        rod_structure.cpp
    )
    target_compile_definitions(ffea_mb PRIVATE FFEA_PARALLEL_PER_BLOB)
    target_link_libraries(ffea_mb PRIVATE ffea_lib)
endif()

#####################################
#### Build Configuration Options ####
#####################################

# This setting enables increased version information printing
set(FFEA_DEVEL OFF)
if(FFEA_DEVEL)
    target_compile_definitions(ffea_lib PUBLIC FFEA_DEVEL)
endif()

##############################
#### LIBRARIES AND CFLAGS ####
##############################

# MPI
if(USE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        # If the MPI installation brings in -flto (i.e. Ubuntu 22.04 libmpich-dev), warn about it and suggest a reconfiguration.
        if(MPI_CXX_COMPILE_OPTIONS MATCHES ".*\-flto.*")
            message(WARNING
                " MPI_CXX_COMPILE_OPTIONS contains '-flto' which is likely to result in linker errors. \n"
                " Consider an alternate MPI implementation which does not embed -flto,\n"
                " Or reconfiguring CMake with -DMPI_CXX_COMPILE_OPTIONS=\"\" if linker error occur.")
        endif()
        # This sets up include directories and lib
        target_link_libraries(ffea_lib PUBLIC MPI::MPI_CXX)
        target_compile_definitions(ffea_lib PUBLIC USE_MPI)
    else()
        message(FATAL_ERROR "MPI C++ Support was not found!")
    endif()
endif()

# Precision
if(USE_PRECISION EQUAL "1")
    target_compile_definitions(ffea_lib PUBLIC USE_DOUBLE_LESS)
elseif(USE_PRECISION EQUAL "2")
    target_compile_definitions(ffea_lib PUBLIC USE_DOUBLE)
elseif(USE_PRECISION EQUAL "3")
    target_compile_definitions(ffea_lib PUBLIC USE_DOUBLE_PLUS)
else()
    message(FATAL_ERROR "USE_PRECISION has to be set to either 1, 2 or 3.")
endif()



# Dynamically create version header
find_package(Git 1.6)
set(VERSION_CMAKECONF "")
if(GIT_FOUND)
    execute_process(
         COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:'%H'
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         OUTPUT_VARIABLE PACKAGE_COMMIT 
         ERROR_VARIABLE PACKAGE_COMMIT_ERROR)

    execute_process(
         COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:'%ad' --date=short
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         OUTPUT_VARIABLE PACKAGE_DATE
         ERROR_VARIABLE PACKAGE_DATE_ERROR)

    execute_process(
         COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         OUTPUT_VARIABLE PACKAGE_BRANCH
         ERROR_VARIABLE PACKAGE_BRANCH_ERROR)

    if(NOT "${PACKAGE_COMMIT}" STREQUAL "${PACKAGE_DATE}")
          target_compile_definitions(ffea_lib PUBLIC USE_CMAKECONF)
         set(VERSION_CMAKECONF "TRUE")
         string(STRIP ${PACKAGE_BRANCH} PACKAGE_BRANCH)
    endif(NOT "${PACKAGE_COMMIT}" STREQUAL "${PACKAGE_DATE}")
endif(GIT_FOUND)

configure_file("${PROJECT_SOURCE_DIR}/include/FFEA_version.h.in"
               "${PROJECT_BINARY_DIR}/include/FFEA_version.h" @ONLY)
include_directories(${PROJECT_BINARY_DIR}/include)




# OpenMP
if(USE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        if(CMAKE_GENERATOR MATCHES "Visual Studio")
            if(MSVC_TOOLSET_VERSION GREATER_EQUAL 142) 
                # Visual Studio 2019+ requires experimental flag for OMP SIMD support
                set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "-openmp:experimental")
                # CMake 3.30 alt: set(OpenMP_RUNTIME_MSVC "experimental")
            else()
                unset(OPENMP_FOUND)
                message(FATAL_ERROR "Visual Studio earlier than 2019 does not support required OpenMP features")
            endif()
        endif()
    else()
        message(FATAL_ERROR "MPI C++ Support was not found!")
    endif()
endif()
if(OPENMP_FOUND)
    target_compile_definitions(ffea_lib PUBLIC USE_OPENMP)
    target_link_libraries(ffea_lib PUBLIC OpenMP::OpenMP_CXX)
endif()


# Parallel futures
# Note, this should now be widely supported can it be made mandatory?
if(USE_FUTURE)
    target_compile_definitions(ffea_lib FFEA_PARALLEL_FUTURE)
endif()


# Math
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(ffea_lib PUBLIC ${MATH_LIBRARY})
endif()


# BOOST
if(USE_BOOST_INTERNAL) 
    message(STATUS "Configuring internal Boost")
    include(../cmake/GetBoostModern.cmake)
else() 
    set(Boost_NO_BOOST_CMAKE ON)
    set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "You could provide valid paths through either:
    -- BOOST_ROOT or
    -- BOOST_INCLUDEDIR and BOOST_LIBRARYDIR or
    try reconfiguring with -DUSE_BOOST_INTERNAL")
endif()
endif()
target_link_libraries(ffea_lib PUBLIC Boost::core Boost::program_options Boost::math Boost::algorithm)

# EIGEN3
if(USE_EIGEN3_INTERNAL) 
    message(STATUS "Configuring internal Eigen")
    # 6.a.1 # set up EIGEN3_HOME
    set(EIGEN3_HOME ${PROJECT_BINARY_DIR}/external/install_eigen)
    # 6.a.2 # get and install Eigen in EIGEN3_HOME
    include(../cmake/GetEigen3Modern.cmake)
else()
    # 6.1 # set up EIGEN3_HOME
    if(EIGEN3_HOME)  
        set(EIGEN3_HOME ${EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
    else()
        set(EIGEN3_HOME $ENV{EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
    endif()
  # 6.2 # we guarantee that if $EIGEN3_HOME is given, it will be take 
        # before than any of the other possible installations available.
  set(CMAKE_PREFIX_PATH ${EIGEN3_HOME} $ENV{EIGEN_HOME} ${CMAKE_PREFIX_PATH})
  # 6.3 # Sometimes "FindEigen3.cmake" is not distributed with CMake, so we're shipping it:
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
  # 6.4 # Go for it:
  find_package(Eigen3 3.2.0 REQUIRED)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"
     AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0"
     AND ${EIGEN3_VERSION} VERSION_LESS "3.3"
     AND NOT ${EIGEN3_VERSION} VERSION_EQUAL "3.2.10")
        message(FATAL_ERROR "\n Found Eigen3 ${EIGEN3_VERSION}, but GCC 5.0 and higher need Eigen version 3.2.10, or > 3.3. Release candidates for 3.3 (internally known as 3.2.91 and higher) are known to fail.")
  endif()
endif() 
target_link_libraries(ffea_lib PUBLIC Eigen3::Eigen)


# Setup install actions
install(TARGETS ffea RUNTIME DESTINATION bin)

if(USE_OPENMP)
    install(TARGETS ffea ffea_mb RUNTIME DESTINATION bin)
endif(USE_OPENMP)

