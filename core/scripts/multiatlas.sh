#!/bin/sh

TARGET=$1
shift

ATLAS_I=""
ATLAS_L=""
while [ "$1" != "" ]; do
    ATLAS_I="${ATLAS_I} $1"
    ATLAS_L="${ATLAS_L} $2"
    shift 2
done

echo ${ATLAS_I}
echo ${ATLAS_L}

tmpdir=`mktemp -d`
echo ${tmpdir}

idx=0
for atlas in ${ATLAS_I}; do
    xfm=${tmpdir}/${idx}.xfm
    cmtk registrationx --echo --init com --auto-multi-levels 5 --omit-original-data --delta-f-threshold 0.01 --dofs 6,9 --ncc --pad-ref 0 --symmetric -o ${xfm} ${TARGET} ${atlas}
    idx=`expr ${idx} + 1`
done

idx=0
for atlas in ${ATLAS_I}; do
    xfm=${tmpdir}/${idx}.xfm
    ffd=${tmpdir}/${idx}.ffd
    cmtk warpx --echo --grid-spacing 80 --grid-refine 3 --min-stepsize 0.25 --max-stepsize 16 --smoothness-constraint-weight 1e-1 --omit-original-data --delta-f-threshold 0.01 --fast --nmi -o ${ffd} --initial ${xfm} ${TARGET} ${atlas}
    idx=`expr ${idx} + 1`
done

idx=0
for atlas in ${ATLAS_I}; do
    ffd=${tmpdir}/${idx}.ffd
    out=${tmpdir}/${idx}_i.nii
    cmtk reformatx --floating ${atlas} --cubic -o ${out} ${TARGET} ${ffd}
    idx=`expr ${idx} + 1`
done

lvote_inputs=""

idx=0
for atlas in ${ATLAS_L}; do
    ffd=${tmpdir}/${idx}.ffd
    out=${tmpdir}/${idx}_l.nii
    cmtk reformatx --floating ${atlas} --pv -o ${out} ${TARGET} ${ffd}

    lvote_inputs="${lvote_inputs} ${tmpdir}/${idx}_i.nii ${tmpdir}/${idx}_l.nii"

    bin/lvote --echo -o output_${idx}_sba.nii --use-sba --patch-radius 5 ${TARGET} ${lvote_inputs}
    
    idx=`expr ${idx} + 1`
done

