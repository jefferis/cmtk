#!/bin/sh

##
##  Copyright 2010-2013 SRI International
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

#
# This script creates a tar/gzip archive of the current CMTK code tree
#

tmpdir=`mktemp -d`

svn export https://nitrc.org/svn/cmtk/trunk/core ${tmpdir}/core

pushd ${tmpdir}
version_major=`fgrep "SET(CMTK_VERSION_MAJOR" core/CMakeLists.txt | sed 's/.*\"\(.*\)\".*/\1/g'`
version_minor=`fgrep "SET(CMTK_VERSION_MINOR" core/CMakeLists.txt | sed 's/.*\"\(.*\)\".*/\1/g'`
version_patch=`fgrep "SET(CMTK_VERSION_PATCH" core/CMakeLists.txt | sed 's/.*\"\(.*\)\".*/\1/g'`

version="${version_major}.${version_minor}.${version_patch}"

eval "gtar -czvf CMTK-${version}-Source.tar.gz --exclude=.svn --transform='s/^core/cmtk-${version}/g' --show-stored-names core"

popd
mv ${tmpdir}/CMTK-${version}-Source.tar.gz .

rm -rf ${tmpdir}
