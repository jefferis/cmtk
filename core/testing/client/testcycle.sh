#!/bin/sh

##
##  Copyright 2009-2010 SRI International
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
##  $Revision: 2602 $
##
##  $LastChangedDate: 2010-12-06 15:23:12 -0800 (Mon, 06 Dec 2010) $
##
##  $LastChangedBy: torstenrohlfing $
##

HOSTNAME=`hostname`
BUILDS=`ls config/${HOSTNAME}`

setup_build()
{
    local build=$1
    if [ ! -d source/${build} ]; then
	svn co https://nitrc.org/svn/cmtk/trunk/ source/${build}
    fi
}

config_build()
{
    local build=$1

    if [ ! -d builds/${build} ]; then
	mkdir -p builds/${build}
    fi

    pushd builds/${build}
    cmake -C ../../config/${HOSTNAME}/${build} ../../source/${build}
    popd
}

CTEST=/usr/bin/ctest

run_ctest()
{
    local build=$1
    shift
    local protocol="$*"

    pushd builds/${build}
    ${CTEST} -D ${protocol}
    popd
}

for b in ${BUILDS}; do
    if [ ! -d builds/${b} ]; then
	setup_build ${b}
	config_build ${b}
	run_ctest $b Experimental
    else
	config_build ${b}
	run_ctest $b Continuous
    fi
done
