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

#-----------------------------------------------------------------------------
# Configure third-party contributions.

OPTION(BUILD_CONTRIB_NCANDA "Build, install, and package script from the National Consortium on Alcohol and Neurodevelopment in Adolescence (N-CANDA)" OFF)

IF(BUILD_CONTRIB_NCANDA)
  ADD_SUBDIRECTORY(ncanda)
ENDIF(BUILD_CONTRIB_NCANDA)


OPTION(BUILD_CONTRIB_MUNGER "Build, install, and package Greg Jefferis' munger script (requires perl)" OFF)

IF(BUILD_CONTRIB_MUNGER)
  ADD_SUBDIRECTORY(munger)
ENDIF(BUILD_CONTRIB_MUNGER)
