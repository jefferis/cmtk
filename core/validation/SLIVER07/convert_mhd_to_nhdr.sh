#!/bin/sh

#
# This script extracts geometry information from a MetaImage (mhd) header
# file and creates a matching NRRD header (nhdr) file.
#

in_mhd=$1
out_nhdr=$2

NDIMS=`fgrep DimSize ${in_mhd} | sed 's/.* = //g'`

MATRIX=`fgrep TransformMatrix ${in_mhd} | sed 's/.* = //g'`
DIR1=`echo ${MATRIX} | gawk '{print $1 "," $2 "," $3}'`
DIR2=`echo ${MATRIX} | gawk '{print $4 "," $5 "," $6}'`
DIR3=`echo ${MATRIX} | gawk '{print $7 "," $8 "," $9}'`
DIRECTIONS="(${DIR1}) (${DIR2}) (${DIR3})"

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
