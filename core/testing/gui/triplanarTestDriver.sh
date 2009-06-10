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

BUILDNAME=$1
BINDIR=$2
DATADIR=$3
RUNTEST=$4

HOSTNAME=`uname -n`
tmpdir=${BINDIR}/../testing/temporary/${HOSTNAME}/${RUNTEST}
mkdir -p ${tmpdir}

BASELINE=${DATADIR}/testing/baseline/${RUNTEST}
if [ -d ${DATADIR}/testing/baseline/${BUILDNAME}/${RUNTEST} ]; then
    BASELINE=${DATADIR}/testing/${BUILDNAME}/${RUNTEST}
fi

cd ${DATADIR}/testing/inputs

run()
{
    cmd=$*

    echo "cd $PWD; $cmd"
    if $cmd; then
	return
    else
	exit 1
    fi
}

run_eval()
{
    cmd=$*

    echo "cd $PWD; $cmd"
    if eval "$cmd"; then
	return
    else
	exit 1
    fi
}

get_unzipped()
{
    if test -f ${1}.gz; then
	tmp=`mktemp`
	gzip -cd ${1}.gz > ${tmp}
	echo ${tmp}
    else
	echo ${1}
    fi
}

check_result()
{
    baseline=`get_unzipped ${BASELINE}/$1`
    result=`get_unzipped ${tmpdir}/$1`

    if test ! -f ${result}; then
	echo "Results file ${result} does not exist"
	exit 1
    fi

    if test ! -f ${baseline}; then
	echo "Baseline file ${baseline} does not exist"
	exit 1
    fi

    echo "diff ${BASELINE}/$1 ${tmpdir}/$1"
    if diff $result $baseline; then
	return
    else
	exit 1
    fi
}

check_results()
{
    for r in $*; do
	check_result $r
    done
}

case ${RUNTEST} in
    TriplanarDefault)
	run ${BINDIR}/triplanar --exec load pat001_pet.hdr export-panel ${tmpdir}/panel.ppm export-axial ${tmpdir}/axial.ppm export-coronal ${tmpdir}/coronal.ppm export-sagittal ${tmpdir}/sagittal.ppm
	for img in panel axial sagittal coronal; do
	    if check_result ${img}.ppm; then
		echo $img ok.
	    else
		exit $?
	    fi
	done
	;;
    TriplanarZoom300)
	run ${BINDIR}/triplanar --exec load pat001_pet.hdr zoom 300 export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarZoom25)
	run ${BINDIR}/triplanar --exec load pat001_pet.hdr zoom 25 export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarSetPixelWindowLevel)
	run ${BINDIR}/triplanar --exec load pat002_ct.hdr crosshair off goto-pixel 32,20,0 window-level 500:0 export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarColormaps)
	run ${BINDIR}/triplanar --exec load parc1.hdr crosshair off window-level 116:58 colormap Labels export-axial ${tmpdir}/labels.ppm colormap Red export-axial ${tmpdir}/red.ppm colormap Rainbow export-axial ${tmpdir}/rainbow.ppm
	for img in labels red rainbow; do
	    if check_result ${img}.ppm; then
		echo $img ok.
	    else
		exit $?
	    fi
	done
	;;
    TriplanarPhantomAx)
	run ${BINDIR}/triplanar --exec load phantom_ax.hdr crosshair off export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarPhantomSa)
	run ${BINDIR}/triplanar --exec load phantom_sa.hdr crosshair off export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarPhantomCo)
	run ${BINDIR}/triplanar --exec load phantom_co.hdr crosshair off export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarPhantomAxNrrd)
	run ${BINDIR}/triplanar --exec load phantom_ax.nhdr crosshair off export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarPhantomSaNrrd)
	run ${BINDIR}/triplanar --exec load phantom_sa.nhdr crosshair off export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
    TriplanarPhantomCoNrrd)
	run ${BINDIR}/triplanar --exec load phantom_co.nhdr crosshair off export-panel ${tmpdir}/panel.ppm
	check_result panel.ppm
	;;
esac

if [ "${tmpdir}" != "" ]; then
    rm -rf ${tmpdir}
fi
