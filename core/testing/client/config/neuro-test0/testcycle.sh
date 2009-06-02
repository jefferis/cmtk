#!/bin/sh

##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
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

export LC_ALL=POSIX

lockfile=${HOME}/testcycle.lock
if test -f ${lockfile}; then
	exit 1
fi
touch ${lockfile}

if [ ! -d ../data ]; then
    pushd ..
    svn co https://www.nitrc.org:443/svn/cmtk/trunk/data/
    popd
else
    svn update ../data
fi

svn update
tests=`ls ctest-*.cmake`

Xpid=""
if [ "${DISPLAY}" == "" ]; then
	Xvfb :1 -screen 0 1024x768x24 -ac &
	Xpid=$!
	export DISPLAY=:1
fi

for t in ${tests}; do
	ctest -S ${t}
done

if [ "${Xpid}" != "" ]; then
  kill ${Xpid}
fi

rm ${lockfile}
