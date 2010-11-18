##
##  Copyright 2009 Torsten Rohlfing
##
##  Copyright 2010 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http://www.nitrc.org/projects/cmtk/
##
##  The Computational Morphometry Toolkit is free software: you can
##  redistribute it and/or modify it under the terms of the GNU General Public
##  License as published by the Free Software Foundation, either version 3 of
##  the License, or (at your option) any later version.
##
##  The Computational Morphometry Toolkit is distributed in the hope that it
##  will be useful, but WITHOUT ANY WARRANTY; without even the implied
##  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with the Computational Morphometry Toolkit.  If not, see
##  <http://www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake build type")

SET(BUILD_APPS "ON" CACHE BOOL "Build command line applications")
SET(BUILD_DOCUMENTATION "ON" CACHE BOOL "Build documentation")

SET(BUILD_GUI "ON" CACHE BOOL "Build GUI components")
SET(BUILD_SHARED_LIBS "OFF" CACHE BOOL "Build shared libraries")
SET(BUILD_TESTING "ON" CACHE BOOL "Build testing components")
SET(BUILD_VALIDATION "ON" CACHE BOOL "Build validation components")

SET(CMTK_BUILD_NRRD "ON" CACHE BOOL "Build with NRRD file format support")
SET(CMTK_BUILD_SMP "ON" CACHE BOOL "Build with SMP parallelism support")
SET(CMTK_BUILD_SQLITE "ON" CACHE BOOL "Build with database support")
SET(CMTK_BUILD_STACKTRACE "ON" CACHE BOOL "Build with stack backtrace printing for crashes")
SET(CMTK_BUILD_UNSTABLE "OFF" CACHE BOOL "Build unstable, experimental code")

SET(CMTK_COORDINATES_DOUBLE "ON" CACHE BOOL "Use double-precision image coordinates")
SET(CMTK_DATA_DOUBLE "ON" CACHE BOOL "Use double-precision data exchange")
SET(CMTK_NUMERICS_DOUBLE "ON" CACHE BOOL "Use double-precision numerics code")
SET(CMTK_SINGLE_COMMAND_BINARY "OFF" CACHE BOOL "Build a single command line binary tool")
SET(CMTK_TESTING_MEMORYCHECK "OFF" CACHE BOOL "Build with support for memory checking")

SET(CMTK_USE_DCMTK "ON" CACHE BOOL "Build with DICOM support")
SET(CMTK_USE_MPI "OFF" CACHE BOOL "Support for Message passing Interface distributed parallelism")
SET(CMTK_USE_OPENMP "ON" CACHE BOOL "Use OpenMP parallelism")
SET(CMTK_USE_QT "ON" CACHE BOOL "Use Gt toolkit for GUI tools")
SET(CMTK_USE_LZMA "OFF" CACHE BOOL "Support for LZMA on-the-fly decompression")
