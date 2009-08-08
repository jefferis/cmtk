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

bindir=@bindir@
registration=${bindir}/registration
reformatx=${bindir}/reformatx
average_images=${bindir}/average_images

tmpdir=`mktemp -d`
if [ "${tmpdir}" == "" ]; then
    echo "ERROR: cannot create tmpdir"
    exit
fi

outfile=$1
reffile=$2
shift 2
fltfiles="$*"

images=${reffile}
for f in ${fltfiles}; do
    name=`basename ${f}`

    echo "Registration: ${name}"
    xform=${tmpdir}/${name}.list
    ${registration} -i -o ${xform} -a 0.01 --dofs 6 -e 2 ${reffile} ${f}

    echo "Reformatting: ${name}"
    nrrd=${tmpdir}/reformat/${name}.nrrd
    ${reformatx} --sinc-cosine --sinc-window-width 5 -o ${nrrd} --floating ${f} ${reffile} ${xform}
    images="${images} ${nrrd}"
done

echo "Averaging"
${average_images} -o ${outfile} ${images}
echo "Done."

##rm -rf ${tmpdir}
