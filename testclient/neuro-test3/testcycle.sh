#!/bin/sh

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

DISABLED_TESTS="OSX-10.5-gcc-Debug OSX-10.4-gcc-Debug"

export LC_ALL=POSIX
export PATH=${PATH}:/opt/local/bin

lockfile=${HOME}/testcycle.lock
if test -f ${lockfile}; then
    exit 1
fi
touch ${lockfile}

svn update
tests=`ls config/*.cmake`
sdks=`ls sdk/*.cmake`

Xpid=""
if [ "${DISPLAY}" == "" ]; then
    /usr/X11/bin/Xvfb :1 -screen 0 1024x768x24 -ac &
    Xpid=$!
    export DISPLAY=:1
fi

for t in ${tests}; do
    for c in ${sdks}; do
	
	if [ ! -d data ]; then
	    svn co https://www.nitrc.org:443/svn/cmtk/trunk/data/
	else
	    svn update data
	fi

	tname=`basename ${c} .cmake`-`basename ${t} .cmake`

	for d in ${DISABLED_TESTS}; do
	    if [ "${tname}" == "${d}" ]; then
		tname=""
	    fi
	done

	if [ "${tname}" != "" ]; then
	    echo "SET(TEST_NAME ${tname})" > /tmp/testfile.cmake
	    cat ${t} ${c} tail.cmake >> /tmp/testfile.cmake
	    
	    /opt/local/bin/ctest -S /tmp/testfile.cmake
	fi
    done
done

if [ "${Xpid}" != "" ]; then
    kill ${Xpid}
fi

rm ${lockfile}
