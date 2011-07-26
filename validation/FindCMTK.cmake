#
# Based on CMake's "FindITK.cmake" script
#

# - Find a CMTK installation or build tree.

# When CMTK is found, the CMTKConfig.cmake file is sourced to setup the
# location and configuration of CMTK.  Please read this file, or
# CMTKConfig.cmake.in from the CMTK source tree for the full list of
# definitions.  Of particular interest is CMTK_USE_FILE, a CMake source file
# that can be included to set the include directories, library directories,
# and preprocessor macros.  In addition to the variables read from
# CMTKConfig.cmake, this find module also defines
#  CMTK_DIR  - The directory containing CMTKConfig.cmake.  
#             This is either the root of the build tree, 
#             or the lib/InsightToolkit directory.  
#             This is the only cache entry.
#   
#  CMTK_FOUND - Whether CMTK was found.  If this is true, 
#              CMTK_DIR is okay.
#
#  USE_CMTK_FILE - The full path to the UseCMTK.cmake file.  
#                 This is provided for backward 
#                 compatability.  Use CMTK_USE_FILE
#                 instead.

#=============================================================================
# Copyright 2001-2010 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2009 Kitware, Inc., Insight Software Consortium
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# ------------------------------------------------------------------------------
# 
# The above copyright and license notice applies to distributions of
# CMake in source and binary form.  Some source files contain additional
# notices of original copyright by their contributors; see each source
# for details.  Third-party software packages supplied with CMake under
# compatible licenses provide their own copyright notices documented in
# corresponding subdirectories.
# 
# ------------------------------------------------------------------------------
# 
# CMake was initially developed by Kitware with the following sponsorship:
# 
#  * National Library of Medicine at the National Institutes of Health
#    as part of the Insight Segmentation and Registration Toolkit (ITK).
# 
#  * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#    Visualization Initiative.
# 
#  * National Alliance for Medical Image Computing (NAMIC) is funded by the
#    National Institutes of Health through the NIH Roadmap for Medical Research,
#    Grant U54 EB005149.
# 
#  * Kitware, Inc.
#=============================================================================

# Use the Config mode of the find_package() command to find CMTKConfig.
# If this succeeds (possibly because CMTK_DIR is already set), the
# command will have already loaded CMTKConfig.cmake and set CMTK_FOUND.
IF(NOT CMTK_FOUND)
  FIND_PACKAGE(CMTK QUIET NO_MODULE
    NAMES CMTK cmtk
    CONFIGS CMTKConfig.cmake
    )
ENDIF()

SET(CMTK_DIR_MESSAGE "Please set CMTK_DIR to the directory containing CMTKConfig.cmake.  This is either the root of the build tree, or PREFIX/lib/cmtk for an installation.")

IF(CMTK_FOUND)
  # Set USE_CMTK_FILE for backward-compatability.
  SET(USE_CMTK_FILE ${CMTK_USE_FILE})
ELSEIF(CMTK_FIND_REQUIRED)
  MESSAGE(FATAL_ERROR ${CMTK_DIR_MESSAGE})
ELSEIF(NOT CMTK_FIND_QUIETLY)
  MESSAGE(STATUS ${CMTK_DIR_MESSAGE})
ENDIF()
