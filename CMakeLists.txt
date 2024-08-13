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

cmake_minimum_required(VERSION 3.18...3.26)
project(FFEA VERSION 2.7.2 LANGUAGES CXX C)

# set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(PACKAGE_NAME "ffea")
#set(VERSION "2") $PROJECT_VERSION_MAJOR
#set(MINOR   "7") $PROJECT_VERSION_MINOR
#set(PATCH   "2") $PROJECT_VERSION_PATCH
#set(PACKAGE_VERSION "${VERSION}.${MINOR}.${PATCH}") $PROJECT_VERSION


#todo check if any of these are now redundant.
option(USE_BOOST_INTERNAL "Use the sources we ship (ON) or your own boost (OFF)" ON)
mark_as_advanced(USE_BOOST_INTERNAL)
option(USE_EIGEN3_INTERNAL "Build Eigen automatically as building FFEA (ON) or use your own install (OFF)" ON)
mark_as_advanced(USE_EIGEN3_INTERNAL)

option(USE_OPENMP "Build with OpenMP support" ON)
option(USE_FAST "Add automatic compiler flags" ON)
mark_as_advanced(USE_FAST)
option(USE_FUTURE "Use C++11 <future> tasking threads" OFF)
mark_as_advanced(USE_FUTURE)
set(USE_PRECISION 2 CACHE STRING "1: Floats + Doubles; 2: Doubles; 3: Doubles + Long Doubles")
mark_as_advanced(USE_PRECISION)

option(USE_MPI "Use MPI parallelisation: EXPERIMENTAL!" OFF)
mark_as_advanced(USE_MPI)
option(USE_PACK "Add generic automatic compiler flags" OFF)
mark_as_advanced(USE_PACK)
set(USE_PACK_DIST "RPM" CACHE STRING "accepted values are RPM and DEB")
mark_as_advanced(USE_PACK_DIST)

# Why is this not an option?
set(BUILD_DOC "TRY" CACHE STRING "ONLY; YES; TRY; NO")
set_property(CACHE BUILD_DOC PROPERTY STRINGS "YES;NO;TRY")

# Docs
include(./cmake/doxygen.cmake)
if (${BUILD_DOC} STREQUAL "ONLY")
  return()
endif ()


#############################################
###### SOURCES TO BUILD AND INSTALL  ########
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
#############################################




# 9 ### PYTHON ######
#find_package(PythonInterp 2.7 REQUIRED)
# python will be used to run some f.
# Now look for NumPy
#execute_process(
#    COMMAND "${PYTHON_EXECUTABLE}" -c "from __future__ import print_function\ntry:\n\t import numpy\n\t print('OK',end=\"\")\nexcept:\n\t print('FAIL',end=\"\")"
#    OUTPUT_VARIABLE NUMPY_FOUND)
#if ("${NUMPY_FOUND}" STREQUAL "FAIL")
#  message(WARNING "NumPy was not found: Several tests won't be run, and the ffeatools script and API will be mostly unusable")
#endif ("${NUMPY_FOUND}" STREQUAL "FAIL")
#if  ("${NUMPY_FOUND}" STREQUAL "OK")
#  message(STATUS "Found NumPy")
#endif  ("${NUMPY_FOUND}" STREQUAL "OK")




# 10 ### FAST! ##########
if (NOT USE_DEBUG)
 if (NOT USE_PACK)
  if (USE_FAST)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.8")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -fipa-pta")
      else(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.8")
        if (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "7.0")
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -fipa-pta -march=native")
        else(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "7.0")
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -march=native")
        endif (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "7.0")
      endif(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.8")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ipo -no-prec-div -fp-model fast=2 -xHost")
    endif()
  endif(USE_FAST)
else (NOT USE_PACK)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -mtune=generic -static -static-libgcc -static-libstdc++")
  else()
      mymessage(FATAL_ERROR "No packaging without GNU compilers")
  endif()
  if (USE_FAST)
      message(STATUS "USE_FAST changed to USE_PACK!")
  endif(USE_FAST)
 endif(NOT USE_PACK)
endif(NOT USE_DEBUG)


# Debug
if (USE_DEBUG)
	# Let the user know what mayhem they may have caused
	message(STATUS "USE_DEBUG has been removed, instead pass -DCMAKE_BUILD_TYPE=Debug at CMake configure time.")
endif (USE_DEBUG)
if(NOT ${GENERATOR_IS_MULTI_CONFIG})
    if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release")
    endif()
endif()
if (CMAKE_BUILD_TYPE EQUAL "Debug")
    target_compile_options(${PROJECT_NAME} $<$<CONFIG:Debug>:-ggdb3>)
endif()

#message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
#get_directory_property( DirDefs COMPILE_DEFINITIONS )
#message(STATUS "Compiler definitions: ${DirDefs}" )

######## END OF LIBRARIES #######
#################################


#######################################
######## INSTALL THE FFEA TOOLS  ######
# 1 - Set up the install python path
include(${CMAKE_ROOT}/Modules/FindPythonInterp.cmake)
set(PYTHONSTUFF "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/ffeatools")
# Debian and derivatives use dist-packages instead of site-packages. 
if (USE_PACK_DIST STREQUAL "DEB")
  set(PYTHONSTUFF "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/dist-packages/ffeatools")
endif(USE_PACK_DIST STREQUAL "DEB")

# 2 - and go through the folders:
add_subdirectory("${PROJECT_SOURCE_DIR}/ffeatools")
#######################################



#####################################################
############ # # # ENABLE TESTS # # # ###############
enable_testing()
# set_tests_properties(ENVIRONMENT OMP_NUM_THREADS=1)
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/consistency")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/physics")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/ffeatools")
add_subdirectory("${PROJECT_SOURCE_DIR}/tests/rods")
set(IGNORE_NUMPY_TESTS "")
if ("${NUMPY_FOUND}" STREQUAL "FAIL")
  set(IGNORE_NUMPY_TESTS 
       python_load_trajectory
       sphere_diffusion_mass sphere_diffusion_mass_check
       sphere_diffusion_nomass sphere_diffusion_nomass_check
       sphere_mass sphere_mass_check 
       sphere_nomass sphere_nomass_check
       cyl_youngs_mod cyl_youngs_mod_check
       cyl_flexrig cyl_flexrig_check)
  message("--------------------------------------------------------")
  message("-----------------WARNING--------------------------------")
  message(STATUS "NumPy could not be found. This means that:")
  message(STATUS "   'ffeatools' - the python package, will fail to import")
  message(STATUS "   'ffeatools' - the script, will only be partially useful")
  message(STATUS "   Several tests have been discarted and will not run.")
  message(STATUS "Your installation will still be complete.")
  message("--------------------------------------------------------")
endif()
configure_file("${PROJECT_SOURCE_DIR}/CTestCustom.in"
                 "${PROJECT_BINARY_DIR}/CTestCustom.cmake" @ONLY)
#####################################################

###############################
##### INSTALL THE LICENSE #####
install(FILES COPYING DESTINATION share/ffea COMPONENT data)

#####################################################
########### # # # CPACK STUFF # # # #################
# FFEA links to Boost statically,
#  so if packagising for different platforms, we only need their gcc version: 
#  gcc 4.8 - Native - Ubuntu 14.04, Mint 17, Centos 7
#          - Works on - Ubuntu 16.04
#  gcc 5.4 - needs Eigen 3.2.10 - 
#          - Ubuntu 16.04 
#  gcc 4.4 - CentOS 6
set(CPACK_PACKAGE_NAME ${PACKAGE_NAME})
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION})
set(CPACK_PACKAGE_ARCHITECTURE "x86_64")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CPACK_PACKAGE_ARCHITECTURE}")
set(CPACK_PACKAGE_VENDOR "University of Leeds, UK.")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FFEA - a toolkit for molecular simulation using Fluctuating Finite Element Method")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_CONTACT "ffea@leeds.ac.uk")

# # # CPACK: DEB or RPM # # # 
if (USE_PACK_DIST STREQUAL "DEB") 
  set(CPACK_GENERATOR "TBZ2" "DEB")
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_amd64")
elseif (USE_PACK_DIST STREQUAL "RPM")
  set(CPACK_GENERATOR "TBZ2" "RPM")
  # set(CPACK_GENERATOR "TXZ" "RPM")
endif()

# # # CPACK: DEB attributes # # # 
# Choose: 
# sudo dpkg -i ffea_1.0_amd64.deb
# sudo apt-get install -f 
# or 
# sudo apt-get install python2.7 python-numpy python-matplotlib libgomp1
# sudo dpkg -i ffea_1.0_amd64.deb
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libgomp1 (>=5.4.0), python2.7, python-numpy, python-matplotlib")
set(CPACK_DEBIAN_PACKAGE_RECOMMENDS "pymol (>=1.8)")

# # # CPACK: RPM attributes # # # 
# sudo yum install numpy python-matplotlib
# sudo rpm -Uvh ffea_1.0_x86_64.rpm
set(CPACK_RPM_PACKAGE_LICENSE "GPLv3")
set(CPACK_RPM_PACKAGE_REQUIRES "numpy, python-matplotlib")
# pymol 1.6 in CentOS 6 is too old for our plugin and CentOS 7 does not provide pymol.


#must come after all cpack settings!
if (USE_PACK)
  include(CPack)
  message(STATUS "Ready for a ${USE_PACK_DIST} package")
endif(USE_PACK)
