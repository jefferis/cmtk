##
##  Copyright 2012 SRI International
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

##
## This file borrows heavily from the analogous InsightToolkit file
##

## FFTW can be compiled and subsequently linked against
## various data types.
## There is a single set of include files, and then muttiple libraries,
## One for each type.  I.e. libfftw.a-->double, libfftwf.a-->float

## The following logic belongs in the individual package
## mark_as_advanced(USE_FFTWD)
## option(USE_FFTWD "Use double precision FFTW if found" ON)
## mark_as_advanced(USE_FFTWF)
## option(USE_FFTWF "Use single precision FFTW if found" ON)

set(FFTW_INC_SEARCHPATH
  /sw/include
  /usr/include
  /usr/local/include
  /usr/include/fftw
  /usr/local/include/fftw
  )

find_path(FFTW_INCLUDE_PATH fftw3.h ${FFTW_INC_SEARCHPATH})

if(FFTW_INCLUDE_PATH)
  set(FFTW_INCLUDE ${FFTW_INCLUDE_PATH})
endif(FFTW_INCLUDE_PATH)

if(FFTW_INCLUDE)
  include_directories( ${FFTW_INCLUDE})
endif(FFTW_INCLUDE)

get_filename_component(FFTW_INSTALL_BASE_PATH ${FFTW_INCLUDE_PATH} PATH)

set(FFTW_LIB_SEARCHPATH
  ${FFTW_INSTALL_BASE_PATH}/lib
  /usr/lib/fftw
  /usr/local/lib/fftw
  )

mark_as_advanced(FFTWD_LIB)
find_library(FFTWD_LIB fftw3 ${FFTW_LIB_SEARCHPATH}) #Double Precision Lib
find_library(FFTWD_THREADS_LIB fftw3_threads ${FFTW_LIB_SEARCHPATH}) #Double Precision Lib only if compiled with threads support

if(FFTWD_LIB)
  set(FFTWD_FOUND 1)
  get_filename_component(FFTW_LIBDIR ${FFTWD_LIB} PATH)
  if(FFTWD_THREADS_LIB)
    set(FFTWD_LIB ${FFTWD_LIB} ${FFTWD_THREADS_LIB} )
  endif(FFTWD_THREADS_LIB)
endif(FFTWD_LIB)

mark_as_advanced(FFTWF_LIB)
find_library(FFTWF_LIB fftw3f ${FFTW_LIB_SEARCHPATH}) #Single Precision Lib
find_library(FFTWF_THREADS_LIB fftw3f_threads ${FFTW_LIB_SEARCHPATH}) #Single Precision Lib only if compiled with threads support

if(FFTWF_LIB)
  set(FFTWF_FOUND 1)
  get_filename_component(FFTW_LIBDIR ${FFTWF_LIB} PATH)
  if(FFTWF_THREADS_LIB)
    set(FFTWF_LIB ${FFTWF_LIB} ${FFTWF_THREADS_LIB} )
  endif(FFTWF_THREADS_LIB)
endif(FFTWF_LIB)
