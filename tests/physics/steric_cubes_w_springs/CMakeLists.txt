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


add_test(NAME steric_n_springs COMMAND ffea cube_x.ffea)
set_tests_properties(steric_n_springs PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

add_executable(test_steric checkIntersections.cpp
               ${PROJECT_SOURCE_DIR}/src/BlobLite.cpp
               )
target_link_libraries(test_steric PRIVATE ffea_lib)

add_test(NAME steric_n_springs_check COMMAND 
         test_steric 2e-9 ../cube_springs_structure/cube_x.node ../cube_springs_structure/cube_x.top 2e-9 ../cube_springs_structure/cube_x.node ../cube_springs_structure/cube_x.top cube_x_trajectory.ftj 10)
set_tests_properties(steric_n_springs_check PROPERTIES DEPENDS "steric_n_springs")




