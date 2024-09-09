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


add_test(NAME squidgy_steric COMMAND ffea sphere_63_120.ffea)
set_tests_properties(squidgy_steric PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

add_executable(test_squidgy checkSquidgySpheres.cpp
               ${PROJECT_SOURCE_DIR}/src/BlobLite.cpp
               )
               
target_compile_features(test_squidgy PUBLIC cxx_std_17)
target_link_libraries(test_squidgy PRIVATE ffea_lib)

add_test(NAME squidgy_steric_check COMMAND 
         test_squidgy 5e-9 ../sphere_63_120_structure/sphere_63_120.node ../sphere_63_120_structure/sphere_63_120.top 5e-9 ../sphere_63_120_structure/sphere_63_120.node ../sphere_63_120_structure/sphere_63_120.top sphere_63_120_trajectory.ftj 332)
set_tests_properties(squidgy_steric_check PROPERTIES DEPENDS "squidgy_steric")




