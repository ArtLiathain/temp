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

cmake_minimum_required(VERSION 2.8.5)
project(FFEA CXX C)
cmake_policy(VERSION 2.8)

# set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(PACKAGE_NAME "ffea")
set(VERSION "2")
set(MINOR   "7")
set(PATCH   "2")
set(PACKAGE_VERSION "${VERSION}.${MINOR}.${PATCH}")
set(FFEA_DEVEL "FALSE")
if (FFEA_DEVEL STREQUAL "TRUE")
  add_definitions(-DFFEA_DEVEL)
endif (FFEA_DEVEL STREQUAL "TRUE")


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

# Debug Build
option(USE_DEBUG "Add debug compiler flags" OFF)

set(BUILD_DOC "TRY" CACHE STRING "ONLY; YES; TRY; NO")

set (LIBS_TO_LINK)

set (NEW_CXX_NEEDED_MSG "You are using ${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}, but need a more recent C++ compiler")
#################################
#### LIBRARIES AND CFLAGS #######

# benchmark ### MPI #########
if (USE_MPI)
  add_definitions(-DUSE_MPI)
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
  include_directories(${MPI_INCLUDE_PATH})
  set(LIBS_TO_LINK ${LIBS_TO_LINK} ${MPI_LIBRARIES}) 
endif(USE_MPI)


include(cmake/Colour.cmake)

# 1 #### PRECISION ###
if (USE_PRECISION EQUAL "1")
  add_definitions(-DUSE_DOUBLE_LESS)
else(USE_PRECISION EQUAL "1")
  if (USE_PRECISION EQUAL "2") 
    add_definitions(-DUSE_DOUBLE)
  else (USE_PRECISION EQUAL "2")
    if (USE_PRECISION EQUAL "3") 
      add_definitions(-DUSE_DOUBLE_PLUS)
    else (USE_PRECISION EQUAL "3")
      mymessage(FATAL_ERROR "USE_PRECISION has to be set to either 1, 2 or 3.")
    endif(USE_PRECISION EQUAL "3")
  endif(USE_PRECISION EQUAL "2")
endif(USE_PRECISION EQUAL "1")
    


# 2 ### CONFIGURE VERSION: 
find_package(Git 1.6)
set(VERSION_CMAKECONF "")
if (GIT_FOUND)
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

  if (NOT "${PACKAGE_COMMIT}" STREQUAL "${PACKAGE_DATE}")
     add_definitions(-DUSE_CMAKECONF)
     set(VERSION_CMAKECONF "TRUE")
     string(STRIP ${PACKAGE_BRANCH} PACKAGE_BRANCH)
  endif (NOT "${PACKAGE_COMMIT}" STREQUAL "${PACKAGE_DATE}")

endif(GIT_FOUND)

configure_file("${PROJECT_SOURCE_DIR}/include/FFEA_version.h.in"
               "${PROJECT_BINARY_DIR}/include/FFEA_version.h" @ONLY)
include_directories(${PROJECT_BINARY_DIR}/include)


# 3 ### DOXYGEN ##### 
if ((${BUILD_DOC} STREQUAL "ONLY") OR (${BUILD_DOC} STREQUAL "YES"))
  message(STATUS "Requested to build the documentation")
  find_package(Doxygen 1.8 REQUIRED)
endif ((${BUILD_DOC} STREQUAL "ONLY") OR (${BUILD_DOC} STREQUAL "YES"))
if (${BUILD_DOC} STREQUAL "TRY")
  find_package(Doxygen 1.8)
endif (${BUILD_DOC} STREQUAL "TRY")

if(DOXYGEN_FOUND)
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/doc/doc/html")
  configure_file("${PROJECT_SOURCE_DIR}/Doxygen.in"
                 "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
  configure_file("${PROJECT_SOURCE_DIR}/DoxygenPyM.in"
                 "${PROJECT_BINARY_DIR}/DoxyfilePyM" @ONLY)
  add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE}
                        ${PROJECT_BINARY_DIR}/DoxyfilePyM
                        DEPENDS docrunner
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/doc/doc/html"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
  add_custom_target(docrunner ALL ${DOXYGEN_EXECUTABLE}
                        ${PROJECT_BINARY_DIR}/Doxyfile
                        WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/doc"
                        COMMENT "Generating documentation for developers"
                        VERBATIM)
  install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/doc/html" DESTINATION
                share/ffea/doc OPTIONAL)
endif(DOXYGEN_FOUND)

if (${BUILD_DOC} STREQUAL "ONLY")
  return()
endif (${BUILD_DOC} STREQUAL "ONLY")



# 4 #### OPENMP #####
if (USE_OPENMP)
  add_definitions(-DUSE_OPENMP)
endif(USE_OPENMP)
include(${CMAKE_ROOT}/Modules/FindOpenMP.cmake)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else(OPENMP_FOUND)
    if (USE_OPENMP)
      mymessage(FATAL_ERROR "Could not find a CXX_FLAG for OpenMP. Possible
solutions:
  -- Disable OPENMP option (-DUSE_OPENMP=0)
  -- Provide valid compiler flags for variable: OpenMP_CXX_FLAGS")
    else (USE_OPENMP)
      mymessage(FATAL_ERROR "Could not find a CXX_FLAG for OpenMP.  Possible
solutions:
  -- Provide valid compiler flags for variable: OpenMP_CXX_FLAGS
  -- Submit a request to the developers suggesting to change omp_get_wtime to something else")
    endif (USE_OPENMP)
endif(OPENMP_FOUND)


# 5 #### USE C++11 ############ 
# g++: -std=c++0x / -std=c++11
# icpc: -std=c++0x / -std=c++11 on LINUX, OS X and WINDOWS 
# pg++: 13.2 crashes with flag -std=c++0x, besides it does not support
#         the latest standard ("cannot open source file "future"). 
# MSVC: does not need any flag.
#include(cmake/CheckForCXX11Features.cmake)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_COMPILER_FLAGS}")
#if (USE_FUTURE)
#  if (HAS_CXX11_FUTURE)
#    message(STATUS "Will use C++11 <future>")
#    add_definitions(-DFFEA_PARALLEL_FUTURE)
#  else (HAS_CXX11_FUTURE)
#    mymessage (FATAL_ERROR "Cannot use C++11 <future>; we did not manage to compile a test program. ${NEW_CXX_NEEDED_MSG} (e. g. GNU 4.8 or higher) or configuring without future: -DUSE_FUTURE=OFF")
#  endif (HAS_CXX11_FUTURE)
#endif(USE_FUTURE)


# 6 ### MATH ######
if (UNIX)
  set(LIBS_TO_LINK ${LIBS_TO_LINK} "m") 
endif(UNIX)


# 7 ### BOOST #####
# set (Boost_NO_SYSTEM_PATHS ON)
if (USE_BOOST_INTERNAL) 

   message(STATUS "Configuring internal Boost")
   if (UNIX) 
      include(cmake/GetBoost.cmake)
   else (UNIX) 
      mymessage(FATAL_ERROR "Some day my prince will come. So far Windows is not supported.")
   endif(UNIX)

else(USE_BOOST_INTERNAL) 
   set (Boost_NO_BOOST_CMAKE ON)
   set (Boost_USE_STATIC_LIBS ON)
   # Check version:
   find_package( Boost 1.54.0 REQUIRED COMPONENTS program_options )
   if (NOT Boost_FOUND)
     mymessage(FATAL_ERROR "You could provide valid paths through either:
      -- BOOST_ROOT or
      -- BOOST_INCLUDEDIR and BOOST_LIBRARYDIR or
try reconfiguring with -DUSE_BOOST_INTERNAL")
   endif (NOT Boost_FOUND)

endif(USE_BOOST_INTERNAL) 

# Anyway:
include_directories(${Boost_INCLUDE_DIRS})
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -L${Boost_LIBRARY_DIRS}")
set(LIBS_TO_LINK ${LIBS_TO_LINK} ${Boost_LIBRARIES})



# 8 ### EIGEN3 #####
if (USE_EIGEN3_INTERNAL) 

   message(STATUS "Configuring internal Eigen")
  # 6.a.1 # set up EIGEN3_HOME
  set(EIGEN3_HOME ${PROJECT_BINARY_DIR}/external/install_eigen)
  # 6.a.2 # get and install Eigen in EIGEN3_HOME
  include(cmake/GetEigen3.cmake)

else (USE_EIGEN3_INTERNAL) 

  # 6.1 # set up EIGEN3_HOME
  if (EIGEN3_HOME)  
    set(EIGEN3_HOME ${EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
  else(EIGEN3_HOME)   
    set(EIGEN3_HOME $ENV{EIGEN3_HOME} CACHE PATH "[OPTIONAL] Path to your local EIGEN3 installation")
  endif(EIGEN3_HOME)

  # 6.2 # we guarantee that if $EIGEN3_HOME is given, it will be take 
        # before than any of the other possible installations available.
  set(CMAKE_PREFIX_PATH ${EIGEN3_HOME} $ENV{EIGEN_HOME} ${CMAKE_PREFIX_PATH})
  # 6.3 # Sometimes "FindEigen3.cmake" is not distributed with CMake, so we're shipping it:
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
  # 6.4 # Go for it:
  find_package(Eigen3 3.2.0 REQUIRED)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
     if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
       if (${EIGEN3_VERSION} VERSION_LESS "3.3") 
         if (NOT ${EIGEN3_VERSION} VERSION_EQUAL "3.2.10")
           mymessage(FATAL_ERROR "\n Found Eigen3 ${EIGEN3_VERSION}, but GCC 5.0 and higher need Eigen version 3.2.10, or > 3.3. Release candidates for 3.3 (internally known as 3.2.91 and higher) are known to fail.")
         endif (NOT ${EIGEN3_VERSION} VERSION_EQUAL "3.2.10")
       endif (${EIGEN3_VERSION} VERSION_LESS "3.3") 
     endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0")
  endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  

endif (USE_EIGEN3_INTERNAL) 
include_directories(${EIGEN3_INCLUDE_DIR})



# 9 ### PYTHON ######
find_package(PythonInterp 2.7 REQUIRED)
# python will be used to run some f.
# Now look for NumPy
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "from __future__ import print_function\ntry:\n\t import numpy\n\t print('OK',end=\"\")\nexcept:\n\t print('FAIL',end=\"\")"
    OUTPUT_VARIABLE NUMPY_FOUND)
if ("${NUMPY_FOUND}" STREQUAL "FAIL")
  message(WARNING "NumPy was not found: Several tests won't be run, and the ffeatools script and API will be mostly unusable")
endif ("${NUMPY_FOUND}" STREQUAL "FAIL")
if  ("${NUMPY_FOUND}" STREQUAL "OK")
  message(STATUS "Found NumPy")
endif  ("${NUMPY_FOUND}" STREQUAL "OK")




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

#
# Debug last, to overwrite everything else
#
if (USE_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)

	# Just in case
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")

	# Let the user know what mayhem they may have caused
	message(STATUS "\n\tDebug mode enabled\n")

endif (USE_DEBUG)

message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
get_directory_property( DirDefs COMPILE_DEFINITIONS )
message(STATUS "Compiler definitions: ${DirDefs}" )

######## END OF LIBRARIES #######
#################################



#############################################
###### SOURCES TO BUILD AND INSTALL  ########
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
#############################################

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
if (${NUMPY_FOUND} STREQUAL "FAIL")
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
endif(${NUMPY_FOUND} STREQUAL "FAIL")
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

