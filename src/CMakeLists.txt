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

add_library(ffea_lib
    FFEA_return_codes.cpp FFEA_input_reader.cpp
    BEM_Poisson_Boltzmann.cpp mesh_node.cpp
    SimulationParams.cpp tetra_element_linear.cpp
    SparseSubstitutionSolver.cpp Face.cpp BiCGSTAB_solver.cpp
    GaussianQuadrature_1d.cpp GaussianQuadrature_tri.cpp
    MassLumpedSolver.cpp MassMatrixLinear.cpp MassMatrixQuadratic.cpp
    MatrixFixedSparsityPattern.cpp LJ_matrix.cpp BindingSite.cpp
    mat_vec_fns_II.cpp NearestNeighbourLinkedListCube.cpp
    PoissonMatrixQuadratic.cpp SparseMatrixTypes.cpp
    SparseMatrixUnknownPattern.cpp SparsityPattern.cpp
    Spring.cpp VdW_solver.cpp Steric_solver.cpp LJSteric_solver.cpp
    GenSoftSSINT_solver.cpp
    KineticState.cpp PreComp_solver.cpp CheckTetrahedraOverlap.cpp
    VolumeIntersection.cpp RngStream.cpp rod_math_v9.cpp
    rod_interactions.cpp rod_structure.cpp rod_blob_interface.cpp
    ffea_test.cpp)

target_compile_features(ffea_lib PUBLIC cxx_std_17)
target_link_libraries(ffea_lib stdc++fs)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ffea ffea.cpp CG_solver.cpp ConjugateGradientSolver.cpp
    FFEA_user_info.cpp mat_vec_fns.cpp NoMassCGSolver.cpp
    SparseMatrixFixedPattern.cpp World.cpp Blob.cpp rod_math_v9.cpp
    rod_interactions.cpp rod_structure.cpp)
if (USE_OPENMP)
  target_compile_definitions(ffea PRIVATE -DFFEA_PARALLEL_WITHIN_BLOB)
endif (USE_OPENMP)

if (USE_OPENMP)
  add_executable(ffea_mb ffea.cpp CG_solver.cpp ConjugateGradientSolver.cpp
    FFEA_user_info.cpp mat_vec_fns.cpp NoMassCGSolver.cpp
    SparseMatrixFixedPattern.cpp World.cpp Blob.cpp  rod_math_v9.cpp
    rod_interactions.cpp rod_structure.cpp)
  target_compile_definitions(ffea_mb PRIVATE -DFFEA_PARALLEL_PER_BLOB)
endif(USE_OPENMP)

add_dependencies(ffea ffea_lib)
if (USE_OPENMP)
  add_dependencies(ffea_mb ffea_lib)
endif(USE_OPENMP)

if (USE_BOOST_INTERNAL)
   add_dependencies(ffea_lib ourboost)
endif(USE_BOOST_INTERNAL)
if (USE_EIGEN3_INTERNAL)
   add_dependencies(ffea_lib getEigen)
endif(USE_EIGEN3_INTERNAL)


target_link_libraries(ffea ffea_lib ${LIBS_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS ffea RUNTIME DESTINATION bin)

if(USE_OPENMP)
    target_link_libraries(ffea_mb ffea_lib ${LIBS_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})
    install (TARGETS ffea ffea_mb RUNTIME DESTINATION bin)
endif(USE_OPENMP)

