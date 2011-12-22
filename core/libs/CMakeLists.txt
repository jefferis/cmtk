##
##  Copyright 1997-2009 Torsten Rohlfing
##
##  Copyright 2004-2011 SRI International
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

SET(BUILD_LIBS_LIST System Numerics Base IO Registration Segmentation Recon)

IF(CMTK_USE_QT)
  SET(BUILD_LIBS_LIST ${BUILD_LIBS_LIST} Pipeline Qt)
ENDIF(CMTK_USE_QT)

IF(CMTK_USE_CUDA)
  SET(BUILD_LIBS_LIST ${BUILD_LIBS_LIST} GPU)
ENDIF(CMTK_USE_CUDA)

IF(CMTK_BUILD_UNSTABLE)
  SET(BUILD_LIBS_LIST ${BUILD_LIBS_LIST} Unstable)
ENDIF(CMTK_BUILD_UNSTABLE)

FOREACH(lib ${BUILD_LIBS_LIST})
  ADD_SUBDIRECTORY(${lib})
ENDFOREACH(lib ${BUILD_LIBS_LIST})
