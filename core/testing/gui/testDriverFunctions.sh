#!/bin/sh

##
##  Copyright 1997-2009 Torsten Rohlfing
##  Copyright 2004-2009 SRI International
##
##  This file is part of the Computational Morphometry Toolkit.
##
##  http:##www.nitrc.org/projects/cmtk/
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
##  <http:##www.gnu.org/licenses/>.
##
##  $Revision$
##
##  $LastChangedDate$
##
##  $LastChangedBy$
##

BUILDNAME=$1
BINDIR=$2
SRCDIR=$3
RUNTEST=$4

HOSTNAME=`uname -n`
tmpdir=${BINDIR}/../testing/temporary/${HOSTNAME}/${RUNTEST}
mkdir -p ${tmpdir}

INPUTS=${SRCDIR}/../../data/testing/

BASELINE=${SRCDIR}/../baseline/${RUNTEST}
if [ -d ${SRCDIR}/../baseline/${BUILDNAME}/${RUNTEST} ]; then
    BASELINE=${SRCDIR}/../baseline/${BUILDNAME}/${RUNTEST}
fi

cd ${INPUTS}

run()
{
    local cmd=$*

    if ! $cmd; then
	exit 1
    fi
}

get_unzipped()
{
    if [ -e ${1}.gz ]; then
	local tmp=`mktemp`
	gzip -cd ${1}.gz > ${tmp}
	echo ${tmp}
    else
	echo ${1}
    fi
}

check_result()
{
    local result=`get_unzipped $1`
    local baseline=`get_unzipped $2`

    echo "diff $1 $2"
    if ! diff $result $baseline; then
	exit 1
    fi
}
