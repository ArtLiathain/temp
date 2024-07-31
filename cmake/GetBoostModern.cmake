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

# Boost+CMake info: https://github.com/boostorg/cmake

set(BOOST_INCLUDE_LIBRARIES "program_options;math;algorithm")
set(BOOST_ENABLE_CMAKE ON)

include(FetchContent)
set(BOOST_DOWNLOAD_VERSION "boost-1.85.0")
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.85.0/boost-1.85.0-cmake.tar.gz
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(Boost)
# Mark some CACHE vars as advanced for a cleaner CMake GUI
mark_as_advanced(FETCHCONTENT_QUIET)
mark_as_advanced(FETCHCONTENT_BASE_DIR)
mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED) 
mark_as_advanced(FETCHCONTENT_SOURCE_DIR_BOOST)
mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_BOOST)

#find_package(Boost REQUIRED COMPONENTS PATHS "${CMAKE_BUILD_DIRECTORY}/_deps")

