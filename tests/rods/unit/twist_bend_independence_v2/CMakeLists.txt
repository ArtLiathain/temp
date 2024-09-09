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

set (TESTROD "${PROJECT_BINARY_DIR}/tests/rods/unit/twist_bend_independence_v2")

add_executable(twist_bend_independence_v2 twist_bend_independence_v2.cpp
               ${PROJECT_SOURCE_DIR}/src/rod_math_v9.cpp # Part of ffea target rather than ffea_lib
               ${PROJECT_SOURCE_DIR}/src/rod_interactions.cpp # Part of ffea target rather than ffea_lib
               ${PROJECT_SOURCE_DIR}/src/rod_structure.cpp # Part of ffea target rather than ffea_lib
               )
target_link_libraries(twist_bend_independence_v2 PRIVATE ffea_lib)

file (COPY twisted_bent_rod.rod DESTINATION ${TESTROD})

add_test(NAME twist_bend_independence_v2 COMMAND twist_bend_independence_v2)
