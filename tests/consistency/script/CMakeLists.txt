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

set (TESTSCRIPT "${PROJECT_BINARY_DIR}/tests/consistency/script")

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(test_script_loader test_script_loader.cpp
               ${PROJECT_SOURCE_DIR}/src/FFEA_user_info.cpp
               )

file (COPY 2r5u_8ang.ffea DESTINATION ${TESTSCRIPT})

target_link_libraries(test_script_loader PUBLIC ffea_lib)

add_test(NAME test_script_loader COMMAND test_script_loader)

