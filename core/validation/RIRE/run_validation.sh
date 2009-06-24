#!/bin/sh

##
##  Copyright 2009 SRI International
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

CMTK_BIN=$1
IN_TREE=$2

PATIENTS=`cd ${IN_TREE} ; ls | fgrep patient`
for p in ${PATIENTS}; do
    mkdir ${p}
    REFSCANS=`cd ${IN_TREE}/${p}; ls | fgrep t`
    MRSCANS=`cd ${IN_TREE}/${p}; ls | fgrep mr`
    for ref in ${REFSCANS}; do	
	for flt in ${MRSCANS}; do
	    ${CMTK_BIN} ${IN_TREE}/${p}/${ref}/header.ascii ${IN_TREE}/${p}/${flt}/header.ascii > ${p}/${p}_${ref}_${flt}
	done
    done
done
