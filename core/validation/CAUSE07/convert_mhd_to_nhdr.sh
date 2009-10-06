#!/bin/sh

#
# This script extracts geometry information from a MetaImage (mhd) header
# file and creates a matching NRRD header (nhdr) file.
#

#
# IMPORTANT:
#
# To get the CAUSE07 images oriented correctly, we need to negate
# the y coordinate direction
#

in_mhd=$1
out_nhdr=$2

NDIMS=`fgrep DimSize ${in_mhd} | sed 's/.* = //g'`

SPACING=`fgrep ElementSpacing ${in_mhd} | sed 's/.* = //g'`
DIR1=`echo ${SPACING} | cut -d" " -f1`
DIR2=`echo ${SPACING} | cut -d" " -f2`
DIR3=`echo ${SPACING} | cut -d" " -f3`
DIRECTIONS="(${DIR1},0,0) (0,-${DIR2},0) (0,0,${DIR3})"

FILE=`fgrep ElementDataFile ${in_mhd} | sed 's/.* = //g'`

echo "NRRD0004
# Complete NRRD file format specification at:
# http://teem.sourceforge.net/nrrd/format.html
type: short
dimension: 3
space: right-anterior-superior
sizes: ${NDIMS}
space directions: ${DIRECTIONS}
kinds: domain domain domain
endian: little
encoding: gzip
space origin: (0,0,0)
data file: ${FILE}.gz" > ${out_nhdr}
