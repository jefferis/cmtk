##
##  Copyright 1997-2010 Torsten Rohlfing
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

set(CTEST_PROJECT_NAME "CMTK")
set(CTEST_NIGHTLY_START_TIME "00:00:00 PST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "neuro.sri.com")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=CMTK")
set(CTEST_DROP_SITE_CDASH TRUE)

SET(CTEST_CUSTOM_WARNING_EXCEPTION
  ${CTEST_CUSTOM_WARNING_EXCEPTION}
  "comparison is always (true|false) due to limited range of data type"
  "warning: iteration variable .* is unsigned"
  "Warning: String literal converted to char. in initialization"
  "sqlite3.c:.*: warning: cast from pointer to integer of different size"
  ".*/Utilities/.*"
)

SET(CTEST_CUSTOM_COVERAGE_EXCLUDE
  ${CTEST_CUSTOM_COVERAGE_EXCLUDE}
  ".*/Utilities/.*"
  )
