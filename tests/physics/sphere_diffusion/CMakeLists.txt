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


#add_test(NAME sphere_diffusion_mass COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_coarse_mass_norestart.ffea)
#set_tests_properties(sphere_diffusion_mass PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
#add_test(NAME sphere_diffusion_mass_check COMMAND ${PYTHON_EXECUTABLE} diffusion_mass.py)
#set_tests_properties(sphere_diffusion_mass_check PROPERTIES DEPENDS sphere_diffusion_mass)

add_test(NAME sphere_diffusion_nomass COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_nomass_norestart.ffea)
set_tests_properties(sphere_diffusion_nomass PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
add_test(NAME sphere_diffusion_nomass_check COMMAND ${PYTHON_EXECUTABLE} diffusion_nomass.py sphere_nomass_norestart.ffea)
set_tests_properties(sphere_diffusion_nomass_check PROPERTIES 
                     DEPENDS sphere_diffusion_nomass
            ENVIRONMENT PYTHONPATH=${PROJECT_SOURCE_DIR}/ffeatools/modules:$ENV{PYTHONPATH})

add_test(NAME sphere_diffusion_mass COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_mass_norestart.ffea)
set_tests_properties(sphere_diffusion_mass PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
add_test(NAME sphere_diffusion_mass_check COMMAND ${PYTHON_EXECUTABLE} diffusion_mass.py sphere_mass_norestart.ffea)
set_tests_properties(sphere_diffusion_mass_check PROPERTIES 
                     DEPENDS sphere_diffusion_mass
            ENVIRONMENT PYTHONPATH=${PROJECT_SOURCE_DIR}/ffeatools/modules:$ENV{PYTHONPATH})




