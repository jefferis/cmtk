#!/bin/sh

##
##  Copyright 2007-2012 SRI International
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
##  $Revision: 4584 $
##
##  $LastChangedDate: 2012-11-12 13:33:41 -0800 (Mon, 12 Nov 2012) $
##
##  $LastChangedBy: torstenrohlfing $
##

##
## Helper functions for file locking and dependency checking using procmail's "lockfile" tool
##

# check whether file needs update due to updated dependencies
CMTK_needs_update()
{
    local target=$1
    
    if [ -e ${target}.lock ]; then
	echo "========================================================================"
	echo "TARGET LOCKED ${target}"
	echo "========================================================================"
	return 1
    fi

    shift
    local sources="$*"
    for source in ${sources}; do
	if [ ! -e ${source} ]; then
	    if [ ! -e ${source}.gz ]; then
		if [ ! -e ${source}.xz ]; then
		    echo "========================================================================"
		    echo "TARGET ${target}"
		    echo "MISSING SOURCE ${source}"
		    echo "========================================================================"
		    
		    return 1
		fi
	    fi
	fi
	
	if [ -e ${source}.lock ]; then
	    echo "========================================================================"
	    echo "TARGET ${target}"
	    echo "SOURCE LOCKED ${source}"
	    echo "========================================================================"
	    
	    return 1
	fi
    done

    if [ ! -e ${target} ]; then
	echo "========================================================================"
	echo "CREATE ${target}"
        echo "========================================================================"

	return 0
    fi

    local source
    for source in ${sources}; do
	if test ! -e ${source}; then
	    if test -e ${source}.gz; then
		source=${source}.gz
	    else
		if test -e ${source}.xz; then
		    source=${source}.xz
		fi
	    fi
	fi
	
	if test ${target} -ot ${source}; then
            echo "========================================================================"
	    echo "UPDATE ${target}"
            echo "DUE TO ${source}"
            echo "========================================================================"

	    return 0
	fi
    done
    
    return 1
}

# NFS-safe (hopefully) file locking.
CMTK_lockfile_create()
{
    local lockfile=$1.lock
    
    mkdir -p `dirname ${lockfile}`
    if lockfile -0 -r0 ${lockfile}; then
	trap "rm -f ${lockfile}; exit" INT TERM EXIT
	return 0
    fi
    
    return 1;
}

CMTK_lockfile_delete()
{
    local file=${1}
    local lockfile=${file}.lock

    local realfile=${file}
    [ -f ${realfile} ] || realfile=${file}.gz
    [ -f ${realfile} ] || realfile=${file}.xz
    
    if [ -f ${realfile} ]; then
	# do not update timestamp if file is older than lock (i.e., existed before and did not get updated at all)
	if test ${lockfile} -ot ${realfile}; then
	    touch --no-create -r ${lockfile} ${realfile}
	fi
    fi
    
    rm -f ${lockfile}
    trap - INT TERM EXIT
}

#
# Combine dependency checking with locking of target
#
CMTK_needs_update_and_lock()
{
    local target=$1
    local realtarget=$target
    [ -e ${realtarget} ] || realtarget=${target}.gz
    [ -e ${realtarget} ] || realtarget=${target}.xz
    [ -e ${realtarget} ] || realtarget=${target}
    
    shift
    if CMTK_needs_update ${realtarget} $*; then
	if CMTK_lockfile_create ${target}; then
	    return 0
	else
	    echo "========================================================================"
	    echo "TARGET LOCKED ${target}"
	    echo "========================================================================"
	fi
    fi
    
    return 1
}
