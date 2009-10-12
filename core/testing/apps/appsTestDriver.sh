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
VALGRIND=$5

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

    echo "pushd $PWD; $cmd; popd"
    if ${VALGRIND} $cmd; then
	return
    else
	exit 1
    fi
}

run_eval()
{
    cmd=$*

    echo "pushd $PWD; $cmd; popd"
    if eval "${VALGRIND} $cmd"; then
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
    AffineRegistrationMrMrMSD)
	run ${BINDIR}/registration -i --dofs 6,9 --msd --match-histograms -o ${tmpdir} pat001_mr_T1.hdr pat002_mr_T2.hdr
	check_result registration
	;;
    ShapeBasedAverage)
	run ${BINDIR}/average_edt -o ${tmpdir}/shape_average.hdr -n 255 --interpolate-image parc1.hdr parc2.hdr parc3.hdr
	check_result shape_average.img
	;;
    ShapeBasedAverageInterpolation)
	run ${BINDIR}/average_edt -o ${tmpdir}/shape_average.hdr -n 255 --interpolate-distance parc12_warp.xform parc13_warp.xform
	check_result shape_average.img
	;;
    ConcatAffineABA)
	run ${BINDIR}/concat_affine -o ${tmpdir}/concat.xform affineA.xform affineB.xform affineA.xform
	check_result concat.xform
	;;
    ConcatAffineABAInvert)
	run ${BINDIR}/concat_affine --invert-output -o ${tmpdir}/concat.xform affineA.xform affineB.xform affineA.xform
	check_result concat.xform
	;;
    ConcatAffineAB1A)
	run ${BINDIR}/concat_affine -o ${tmpdir}/concat.xform affineA.xform --inverse affineB.xform affineA.xform
	check_result concat.xform
	;;
    ConcatAffineAA1)
	run ${BINDIR}/concat_affine -o ${tmpdir}/concat.xform affineA.xform --inverse affineA.xform
	check_result concat.xform
	;;
    ConcatAffineA1A)
	run ${BINDIR}/concat_affine -o ${tmpdir}/concat.xform -- --inverse affineA.xform affineA.xform
	check_result concat.xform
	;;
    CongealFromInit)
	run ${BINDIR}/congeal -v --force-background 0 -O ${tmpdir} --dofs 6 --dofs 9 -e 4 -a 0.25 --downsample-from 2 --downsample-to 1 --sampling-density 1.0 --zero-sum groupwise_init_brain123.xforms
	check_result congeal.xforms
	check_result average_congeal.img
	;;
    CongealFromInitSampling)
	run ${BINDIR}/congeal -v --force-background 0 -O ${tmpdir} --dofs 6 --dofs 9 -e 4 -a 0.25 --downsample-to 1 --sampling-density 0.5 --zero-sum groupwise_init_brain123.xforms
	## no baseline; this mode uses probabilistic sampling
	;;
    CongealBackground)
	run ${BINDIR}/congeal -v --force-background 0 -O ${tmpdir} --template spgr_brain_1.hdr --dofs 6 --dofs 9 -e 2 -a 0.25 --downsample-from 2 --downsample-to 1 --sampling-density 1.0 --zero-sum spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_result congeal.xforms
	check_result average_congeal.img
	;;
    CongealUseTemplate)
	run ${BINDIR}/congeal -v --force-background 0 -O ${tmpdir} --template spgr_brain_1.hdr --dofs 6 --dofs 9 -e 2 -a 0.25 --downsample-from 2 --downsample-to 1 --sampling-density 1.0 --template-with-data spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_result congeal.xforms
	check_result average_congeal.img
	;;
    CongealZeroSumSmooth)
	run ${BINDIR}/congeal -v --smooth 0.5 --downsample-to 0 -O ${tmpdir} --template spgr_brain_1.hdr --dofs 6 --dofs 9 -e 2 -a 0.25 --downsample-from 2 --sampling-density 1.0 --zero-sum spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_result congeal.xforms
	check_result average_congeal.img
	;;
    CongealWarpFromInit)
	run ${BINDIR}/congeal_warp -v --force-background 0 -O ${tmpdir} --grid-spacing 200 --refine-grid 1 -e 2 -a 1 --downsample-from 2 --downsample-to 1 groupwise_init_brain123.xforms
	check_results congeal_warp.xforms average_congeal_warp.img
	;;
    CongealWarpFitFromInit)
	run ${BINDIR}/congeal_warp -v --force-background 0 -O ${tmpdir} --grid-spacing 200 --grid-spacing-fit --refine-grid 1 -e 4 -a 1 --downsample-from 2 --downsample-to 1 groupwise_init_brain123.xforms
	check_results congeal_warp.xforms average_congeal_warp.img
	;;
    CongealWarpUseTemplate)
	run ${BINDIR}/congeal_warp -v --force-background 0 -O ${tmpdir} --grid-spacing 200 --refine-grid 1 -e 2 -a 1 --downsample-from 2 --downsample-to 2 --template-with-data spgr_brain_1.hdr groupwise_init_brain123.xforms
	check_results congeal_warp.xforms average_congeal_warp.img
	;;
    CongealWarpFromInitZeroSum)
	run ${BINDIR}/congeal_warp -v --zero-sum --force-background 0 -O ${tmpdir} --grid-spacing 200 --refine-grid 1 -e 2 -a 1 --downsample-from 2 --downsample-to 1 groupwise_init_brain123.xforms
	check_results congeal_warp.xforms average_congeal_warp.img
	;;
    ConvertAutoCrop)
	run_eval "${BINDIR}/convert -v --crop-xform-out ${tmpdir}/crop.xform --auto-crop 1 spgr_brain_1.hdr ${tmpdir}/cropped.nii  > ${tmpdir}/crop.txt"
	cat ${tmpdir}/crop.xform/registration | sed 's/reference_study ".*"/reference_study ""/g' > ${tmpdir}/crop.xform/registration.tmp && mv ${tmpdir}/crop.xform/registration.tmp ${tmpdir}/crop.xform/registration
	check_results cropped.nii crop.txt crop.xform/registration
	;;
    ConvertBoundaryMap)
	run ${BINDIR}/convert --bmap parc1.hdr ${tmpdir}/boundary_map.hdr
	check_result boundary_map.img
	;;
    ConvertBoundaryMapMultiValue)
	run ${BINDIR}/convert --bmap-multi parc1.hdr ${tmpdir}/boundary_map_multi.hdr
	check_result boundary_map_multi.img
	;;
    ConvertDownsample)
	run ${BINDIR}/convert --downsample 8,4,1 phantom_ax.hdr ${tmpdir}/downsample.hdr
	check_results downsample.hdr downsample.img
	;;
    ConvertDownsampleNrrd)
	run ${BINDIR}/convert --downsample 8,4,1 phantom_ax.nhdr ${tmpdir}/downsample.nhdr
	check_results downsample.nhdr downsample.raw
	;;
    ConvertErodeDilateErode)
	run ${BINDIR}/convert --erode 1 --dilate 2 --erode 1 parc1_bin.hdr ${tmpdir}/parc1_ede.img
	check_results parc1_ede.img
	;;
    ConvertDilateErodeDilate)
	run ${BINDIR}/convert --dilate 1 --erode 2 --dilate 1 parc1_bin.hdr ${tmpdir}/parc1_ded.img
	check_results parc1_ded.img
	;;
    ConvertMedianFilter1)
	run ${BINDIR}/convert --median 1 spgr_brain_1.hdr ${tmpdir}/spgr_brain_1.hdr
	check_result spgr_brain_1.img
	;;
    ConvertMedianFilter2)
	run ${BINDIR}/convert --median 2 spgr_brain_1.hdr ${tmpdir}/spgr_brain_1.hdr
	check_result spgr_brain_1.img
	;;
    ConvertNiftiToAnalyze)
	run ${BINDIR}/convert spgr_brain_1.nii ${tmpdir}/spgr_brain_1.hdr
	check_results spgr_brain_1.hdr spgr_brain_1.img
	;;
    ConvertAnalyzeToNifti)
	run ${BINDIR}/convert spgr_brain_1.hdr ${tmpdir}/spgr_brain_1.nii
	check_result spgr_brain_1.nii
	;;
    ConvertNiftiDetachedToNifti)
	run ${BINDIR}/convert spgr_brain_nifti.hdr ${tmpdir}/spgr_brain_1.nii
	check_result spgr_brain_1.nii
	;;
    ConvertAnalyzeToNiftiDetached)
	run ${BINDIR}/convert spgr_brain_1.hdr ${tmpdir}/spgr_brain_1.img
	check_results spgr_brain_1.hdr spgr_brain_1.img
	;;
    Dcm2Image)
	run ${BINDIR}/dcm2image -O ${tmpdir}/%03d.hdr ${PWD}/dcm/
	check_results 001.hdr 001.img 002.hdr 002.img 003.hdr 003.img
	;;
    Dcm2ImageZ)
	run ${BINDIR}/dcm2image -O ${tmpdir}/%03d.nii ${PWD}/dcmz/
	check_results 001.nii 002.nii 003.nii
	;;
    Dcm2ImageNrrd)
	run ${BINDIR}/dcm2image -O ${tmpdir}/%03d.nhdr ${PWD}/dcm/
	check_results 001.nhdr 001.raw 002.nhdr 002.raw 003.nhdr 003.raw
	;;
    DescribeMR1)
	run_eval "${BINDIR}/describe -m parc1_bin.hdr > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    DescribeMR2)
	run_eval "${BINDIR}/describe -m phantom_ax.hdr phantom_co.hdr phantom_sa.hdr > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    DescribeMR3)
	run_eval "${BINDIR}/describe -m header_only.hdr > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    DescribeMR4)
	run_eval "${BINDIR}/describe -m phantom_ax.nii phantom_co.nii phantom_sa.nii > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    DescribeMRBiorad)
	run_eval "${BINDIR}/describe -m bioradvol.PIC.gz > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    DescribeMRNrrd1)
	run_eval "${BINDIR}/describe -m phantom_ax.nhdr phantom_co.nhdr phantom_sa.nhdr > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    DescribeMRNrrd2)
	run_eval "${BINDIR}/describe -m split_ax_0.nhdr split_ax_1.nhdr split_ax_2.nhdr > ${tmpdir}/describe.txt"
	check_result describe.txt
	;;
    FilterGaussian)
	run ${BINDIR}/filter --gaussian 1 --radius 2 rat_fse_erly.hdr ${tmpdir}/filter.hdr
	check_result filter.img
	;;
    FilterGaussianSmallKernel)
	run ${BINDIR}/filter --gaussian 1 --radius 1.1 rat_fse_erly.hdr ${tmpdir}/filter.hdr
	check_result filter.img
	;;
    FilterGaussianNoFilter)
	run ${BINDIR}/filter --gaussian 1 --radius 0.5 rat_fse_erly.hdr ${tmpdir}/filter.hdr
	check_result filter.img
	;;
    FilmFourthOrder)
	run ${BINDIR}/film --coronal --nmi --fourth-order-error --passes 3 --num-iterations 3 --injection-kernel-radius 2 --injection-kernel-sigma 1 --write-injected-image ${tmpdir}/injected.hdr interleave_thnfse.hdr ${tmpdir}/corrected.hdr
	check_result corrected.img
	;;
    FilmCubic)
	run ${BINDIR}/film --coronal --nmi --cubic --passes 3 --num-iterations 3 --injection-kernel-radius 2 --injection-kernel-sigma 1 --write-injected-image ${tmpdir}/injected.hdr interleave_thnfse.hdr ${tmpdir}/corrected.hdr
	check_result corrected.img
	;;
    FilmMSDLinearNoTrunc)
	run ${BINDIR}/film --msd --coronal --linear --no-truncation --passes 3 --num-iterations 3 --injection-kernel-radius 2 --injection-kernel-sigma 1 interleave_thnfse.hdr ${tmpdir}/corrected.hdr
	check_result corrected.img
	;;
    FilmMISincRefSPGR)
	run ${BINDIR}/film --mi --coronal --cubic --passes 3 --num-iterations 3 --injection-kernel-radius 2 --injection-kernel-sigma 1 --reference-image interleave_3dspgr.hdr interleave_thnfse.hdr ${tmpdir}/corrected.hdr
	check_result corrected.img
	;;
    GregxformFordwardBackward)
	run_eval "cat vol001_t0_points.xyz | ${BINDIR}/gregxform -f vol001_mr_t0t1_warp.xform | ${BINDIR}/gregxform vol001_mr_t0t1_warp.xform > ${tmpdir}/vol001_t0_points.xyz"
	check_result vol001_t0_points.xyz
	;;
    GregxformAffine)
	run_eval "cat vol001_t0_points.xyz | ${BINDIR}/gregxform -f vol001_mr_t0t1.list > ${tmpdir}/vol001_t0_points.xyz"
	check_result vol001_t0_points.xyz
	;;
    GregxformAffineFromWarp)
	run_eval "cat vol001_t0_points.xyz | ${BINDIR}/gregxform --affine -f vol001_mr_t0t1_warp.xform > ${tmpdir}/vol001_t0_points.xyz"
	check_result vol001_t0_points.xyz
	;;
    GregxformAffineFromWarpFwdBwd)
	run_eval "cat vol001_t0_points.xyz | ${BINDIR}/gregxform --affine -f vol001_mr_t0t1_warp.xform | ${BINDIR}/gregxform --affine vol001_mr_t0t1_warp.xform > ${tmpdir}/vol001_t0_points.xyz"
	check_result vol001_t0_points.xyz
	;;
    GroupwiseInitCentersOfMass)
	run ${BINDIR}/groupwise_init -O ${tmpdir} --align-centers-of-mass spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_result groupwise_init.xforms
	# cannot compare separate xforms due to local file system path in the file
	check_result groupwise_init_average.img
	;;
    GroupwiseInitCentersOfMassTemplate)
	run ${BINDIR}/groupwise_init -O ${tmpdir} --template spgr_brain_1.hdr --align-centers-of-mass spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_result groupwise_init.xforms
	check_results groupwise_init_pairs/target-000.list/registration groupwise_init_pairs/target-001.list/registration groupwise_init_pairs/target-002.list/registration
	check_result groupwise_init_average.img
	;;
    GroupwiseRMIFromInit)
	run ${BINDIR}/groupwise_rmi --force-background 0 -O ${tmpdir} --dofs 6 --dofs 9 -e 4 -a 0.25 --downsample-from 2 --downsample-to 1 --sampling-density 1.0 --zero-sum groupwise_init_brain123.xforms
	check_results groupwise_rmi.xforms average_groupwise_rmi.img
	;;
    GroupwiseRMIFromInitDeltaF)
	run ${BINDIR}/groupwise_rmi --force-background 0 -O ${tmpdir} --dofs 6,9 -e 4 --delta-f-threshold 0.1 -a 0.25 --downsample-from 2 --downsample-to 1 --sampling-density 1.0 --zero-sum groupwise_init_brain123.xforms
	check_results groupwise_rmi.xforms average_groupwise_rmi.img
	;;
    GroupwiseRMIFromInitSampling)
	run ${BINDIR}/groupwise_rmi --force-background 0 -O ${tmpdir} --dofs 6 --dofs 9 -e 4 -a 0.25 --downsample-to 1 --sampling-density 0.5 --zero-sum groupwise_init_brain123.xforms
	## no baseline; this mode uses probabilistic sampling
	;;
    GroupwiseRMIBackground)
	run ${BINDIR}/groupwise_rmi --force-background 0 -O ${tmpdir} --template spgr_brain_1.hdr --dofs 6 --dofs 9 -e 2 -a 0.25 --downsample-from 2 --downsample-to 1 --sampling-density 1.0 --zero-sum spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_results groupwise_rmi.xforms average_groupwise_rmi.img
	;;
    GroupwiseRMIZeroSumSmooth)
	run ${BINDIR}/groupwise_rmi --smooth 0.5  --downsample-to 0 -O ${tmpdir} --template spgr_brain_1.hdr --dofs 6 --dofs 9 -e 2 -a 0.25 --downsample-from 2 --sampling-density 1.0 --zero-sum spgr_brain_1.hdr spgr_brain_2.hdr spgr_brain_3.hdr
	check_result groupwise_rmi.xforms
	check_result average_groupwise_rmi.img
	;;
    Histogram)
	run ${BINDIR}/histogram -o ${tmpdir}/histogram spgr_3t.hdr
	check_result histogram
	;;
    HistogramNorm)
	run ${BINDIR}/histogram --normalize -o ${tmpdir}/histogram spgr_3t.hdr
	check_result histogram
	;;
    HistogramBinsMinMax)
	run ${BINDIR}/histogram --min 300 --max 3000 --nbins 16 -o ${tmpdir}/histogram spgr_3t.hdr
	check_result histogram
	;;
    HistogramBinsMinMaxTrunc)
	run ${BINDIR}/histogram --truncate --min 300 --max 3000 --nbins 16 -o ${tmpdir}/histogram spgr_3t.hdr
	check_result histogram
	;;
    HistogramMask)
	run ${BINDIR}/histogram --mask spgr_3t_mask.hdr -o ${tmpdir}/histogram spgr_3t.hdr
	check_result histogram
	;;
    ImagemathAverage)
	run ${BINDIR}/imagemath --in parc1_bin.hdr parc2_bin.hdr parc3_bin.hdr --average --out ${tmpdir}/average.hdr
	check_result average.img
	;;
    ImagemathXor)
	run ${BINDIR}/imagemath --in parc1.hdr --scalar-xor 1 --out ${tmpdir}/xor.hdr
	check_result xor.img
	;;
    ImagemathContractLabels)
	run ${BINDIR}/imagemath --in parc1.hdr parc2.hdr parc3.hdr --contract-labels --out ${tmpdir}/contract.hdr
	check_result contract.img
	;;
    ImagemathVote)
	run ${BINDIR}/imagemath --in parc1_bin.hdr parc2_bin.hdr parc3_bin.hdr --vote --out ${tmpdir}/vote.hdr
	check_result vote.img
	;;
    ImagemathStackEntropyLabels)
	run ${BINDIR}/imagemath --in parc1_bin.hdr parc2_bin.hdr parc3_bin.hdr --stack-entropy-labels --out ${tmpdir}/entropy.hdr
	check_result entropy.img
	;;
    ImagemathSTAPLE)
	run ${BINDIR}/imagemath --in parc1_bin.hdr parc2_bin.hdr parc3_bin.hdr --staple 10 --out ${tmpdir}/staple.hdr
	check_result staple.img
	;;
    ImagemathMultiClassSTAPLE)
	run ${BINDIR}/imagemath --in parc1.hdr parc2.hdr parc3.hdr --mstaple 2 --out ${tmpdir}/mstaple.hdr
	check_result mstaple.img
	;;
    ImagemathCombinePCA)
	run ${BINDIR}/imagemath --in rat_fse_erly.hdr rat_fse_late.hdr --combine-pca --trunc --out ${tmpdir}/combine_pca.hdr
	check_result combine_pca.img
	;;
    ImagemathT2)
	run ${BINDIR}/imagemath --float --in rat_fse_erly.hdr --log --in rat_fse_late.hdr --log --scalar-mul -1 --add --one-over --out ${tmpdir}/t2.hdr
	check_result t2.img
	;;
    ImagemathLogOddsAdd)
	run ${BINDIR}/imagemath --float --in pbmap_wm_2.nii --logit --in pbmap_wm_1.nii --logit --average --logistic --out ${tmpdir}/logodds_add.hdr
	check_result logodds_add.img
	;;
    ImagemathLogOddsAdd2)
	run ${BINDIR}/imagemath --float --in pbmap_wm_2.nii pbmap_wm_1.nii --logit-all --average --logistic-all --out ${tmpdir}/logodds_add.hdr
	check_result logodds_add.img
	;;
    ImagemathMatchHistograms)
        run ${BINDIR}/imagemath --in spgr_brain_{1,2}.hdr --match-histograms --out ${tmpdir}/match_histograms.hdr
        check_result match_histograms.img
        ;;
    ImagemathMatchHistogramsPadding)
        run ${BINDIR}/imagemath --set-padding-value 0 --in spgr_brain_{1,2}.hdr --match-histograms --out ${tmpdir}/match_histograms.hdr
        check_result match_histograms.img
        ;;
    LevelsetDefault)
	run ${BINDIR}/levelset -v vol001_mr_t1.hdr ${tmpdir}/levelset.hdr
	check_result levelset.img
	;;
    MakeInitialAffineCenterOfMass)
	run ${BINDIR}/make_initial_affine --xform-ras --mode centers-of-mass box1.hdr box3.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffinePrincipalAxes1)
	run ${BINDIR}/make_initial_affine --xform-ras --mode principal-axes box1.hdr box2.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffinePrincipalAxes2)
	run ${BINDIR}/make_initial_affine --xform-ras --mode principal-axes box1.hdr box3.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffinePrincipalAxes3)
	run ${BINDIR}/make_initial_affine --xform-ras --mode principal-axes box2.hdr box3.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffinePrincipalAxes4)
	run ${BINDIR}/make_initial_affine --xform-ras --mode principal-axes box1.hdr box4.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffinePrincipalAxes5)
	run ${BINDIR}/make_initial_affine --xform-ras --mode principal-axes box2.hdr box4.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffinePrincipalAxes6)
	run ${BINDIR}/make_initial_affine --xform-ras --mode principal-axes box3.hdr box4.hdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffineDirectionVectorsNrrdAxSa)
	run ${BINDIR}/make_initial_affine --xform-ras phantom_ax.nhdr phantom_sa.nhdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffineDirectionVectorsNrrdAxCo)
	run ${BINDIR}/make_initial_affine --xform-ras phantom_ax.nhdr phantom_co.nhdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffineDirectionVectorsNrrdSaCo)
	run ${BINDIR}/make_initial_affine --xform-ras phantom_sa.nhdr phantom_co.nhdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffineDirectionVectorsNrrdAxSaNative)
	run ${BINDIR}/make_initial_affine phantom_ax.nhdr phantom_sa.nhdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffineDirectionVectorsNrrdAxCoNative)
	run ${BINDIR}/make_initial_affine phantom_ax.nhdr phantom_co.nhdr ${tmpdir}/xform
	check_result xform
	;;
    MakeInitialAffineDirectionVectorsNrrdSaCoNative)
	run ${BINDIR}/make_initial_affine phantom_sa.nhdr phantom_co.nhdr ${tmpdir}/xform
	check_result xform
	;;
    McAffine1)
	run ${BINDIR}/mcaffine --downsample-from 4 --downsample-to 1 --initial-step-size 1 --final-step-size 0.5 --dofs 6 --covariance -o ${tmpdir}/xform pat001_mr_T1.hdr -- pat001_pet.hdr
	check_result xform
	;;
    McAffine2)
	run ${BINDIR}/mcaffine --initial-xform McAffine2_initial.xform --downsample-from 4 --downsample-to 1 --initial-step-size 1 --final-step-size 0.5 --dofs 6 --histograms -o ${tmpdir}/xform pat001_mr_T1.hdr -- pat001_pet.hdr
	check_result xform
	;;
    McAffine3)
	run ${BINDIR}/mcaffine --downsample-from 4 --downsample-to 1 --initial-step-size 1 --final-step-size 0.5 --dofs 6 --dofs 9 --covariance -o ${tmpdir}/xform rat_fse_erly.hdr rat_fse_late.hdr -- rat2_fse_erly.hdr rat2_fse_late.hdr
	check_result xform
	;;
    McAffine4)
	run ${BINDIR}/mcaffine --downsample-from 4 --downsample-to 1 --delta-f-threshold 0.1 --initial-step-size 1 --final-step-size 0.5 --dofs 6,9 --covariance -o ${tmpdir}/xform rat_fse_erly.hdr rat_fse_late.hdr -- rat2_fse_erly.hdr rat2_fse_late.hdr
	check_result xform
	;;
    McWarp1)
	# Test breaks when using multiple threads due to floating point effects
	export CMTK_NUM_THREADS=1
	run ${BINDIR}/mcwarp --downsample-from 2 --downsample-to 1 --initial-step-size 1 --final-step-size 0.5 --grid-spacing 14 --refine-grid 2 --covariance -o ${tmpdir}/xform McAffine_rat_rat2.xform
	check_result xform
	unset CMTK_NUM_THREADS
	;;
    McWarp2)
	run ${BINDIR}/mcwarp --downsample-from 2 --downsample-to 1 --initial-step-size 1 --final-step-size 0.5 --grid-spacing 14 --refine-grid 2 --histograms -o ${tmpdir}/xform McAffine_rat_rat2.xform
	check_result xform
	;;
    McWarp3)
	# Test breaks when using multiple threads due to floating point effects
	export CMTK_NUM_THREADS=1
	run ${BINDIR}/mcwarp --downsample-from 2 --downsample-to 1 --delta-f-threshold 0.1 --initial-step-size 1 --final-step-size 0.5 --grid-spacing 14 --refine-grid 2 --covariance -o ${tmpdir}/xform McAffine_rat_rat2.xform
	check_result xform
	unset CMTK_NUM_THREADS
	;;
    MkPhantomBox)
	run ${BINDIR}/mk_phantom_3d -o ${tmpdir}/phantom.nii --dims 10,10,10 --voxel 1,1,1 box 2,2,2 5,5,5 10
	check_results phantom.nii
	;;
    MkPhantomSphere)
	run ${BINDIR}/mk_phantom_3d -o ${tmpdir}/phantom.nii --dims 10,10,10 --voxel 1,1,1 sphere 7,7,7 5 20
	check_results phantom.nii
	;;
    MkPhantomBoxSphere)
	run ${BINDIR}/mk_phantom_3d -o ${tmpdir}/phantom.nii --char --bg 50 --dims 10,10,10 --voxel 1,1,1 sphere 7,7,7 5 20 box 2,2,2 5,5,5 10
	check_results phantom.nii
	;;
    MkPhantomImport)	
	run ${BINDIR}/mk_phantom_3d -o ${tmpdir}/phantom.nii --import spgr_3t_mask.hdr sphere 30,40,30 20 0
	check_results phantom.nii
	;;
    GLMDefault)
	run_eval "${BINDIR}/glm -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img  model_tstat_00_age.img model_param_00_age.img
	check_results model_tstat_01_sex.img model_param_01_sex.img model_tstat_02_CONST.img model_param_02_CONST.img
	;;
    GLMNormalize)
	run_eval "${BINDIR}/glm --normalize -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img  model_tstat_00_age.img model_param_00_age.img
	check_results model_tstat_01_sex.img model_param_01_sex.img model_tstat_02_CONST.img model_param_02_CONST.img
	;;
    GLMExp)
	run_eval "${BINDIR}/glm --exp -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img  model_tstat_00_age.img model_param_00_age.img
	check_results model_tstat_01_sex.img model_param_01_sex.img model_tstat_02_CONST.img model_param_02_CONST.img
	;;
    GLMNoConstant)
	run_eval "${BINDIR}/glm --exclude-constant -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img  model_tstat_00_age.img model_param_00_age.img
	check_results model_tstat_01_sex.img model_param_01_sex.img
	;;
    GLMIgnore)
	run_eval "${BINDIR}/glm --ignore-parameter 1 -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img  model_tstat_00_age.img model_param_00_age.img
	check_results model_tstat_02_CONST.img model_param_02_CONST.img
	;;
    GLMSelect)
	run_eval "${BINDIR}/glm --select-parameter age -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img
	check_results model_tstat_00_age.img model_param_00_age.img model_tstat_02_CONST.img model_param_02_CONST.img
	;; 
    GLMCrop)
	run_eval "${BINDIR}/glm --crop 10,10,0,20,20,2 -v -O ${tmpdir}/model_%s_%02d_%s.hdr jacobian_glm.txt jacobian-%s.nii > ${tmpdir}/model_stdout.txt"
	check_results model_stdout.txt model_fstat_00_model.img  model_tstat_00_age.img model_param_00_age.img
	check_results model_tstat_01_sex.img model_param_01_sex.img model_tstat_02_CONST.img model_param_02_CONST.img
	;;
    MrbiasMulIncremental)
	run ${BINDIR}/mrbias --incremental -M 2 --write-bias-mul ${tmpdir}/bias_mul.hdr --thresh-min 100 spgr_3t.hdr ${tmpdir}/corrected.hdr
	check_results corrected.img bias_mul.img
	;;
    MrbiasMulLogIntensity)
	run ${BINDIR}/mrbias --log-intensities -M 2 --write-bias-mul ${tmpdir}/bias_mul.hdr --thresh-min 100 spgr_3t.hdr ${tmpdir}/corrected.hdr
	check_results corrected.img bias_mul.img
	;;
    MrbiasAddMulMask)
	export CMTK_NUM_THREADS=1
	run ${BINDIR}/mrbias -A 1 -M 2 --mask spgr_3t_mask.hdr --write-bias-add ${tmpdir}/bias_add.hdr --write-bias-mul ${tmpdir}/bias_mul.hdr spgr_3t.hdr ${tmpdir}/corrected.hdr
	unset CMTK_NUM_THREADS
	check_results corrected.img bias_mul.img bias_add.img
	;;
    Overlap)
	run_eval "${BINDIR}/overlap parc1.hdr parc2.hdr parc3.hdr > ${tmpdir}/overlap.txt"
	check_result overlap.txt
	;;
    OverlapNumLabels)
	run_eval "${BINDIR}/overlap -N 2 parc1.hdr parc2.hdr parc3.hdr > ${tmpdir}/overlap.txt"
	check_result overlap.txt
	;;
    OverlapByLabel)
	run_eval "${BINDIR}/overlap --by-label parc1.hdr parc2.hdr parc3.hdr > ${tmpdir}/overlap.txt"
	check_result overlap.txt
	;;
    OverlapFirst)
	run_eval "${BINDIR}/overlap --first-label 10 parc1.hdr parc2.hdr parc3.hdr > ${tmpdir}/overlap.txt"
	check_result overlap.txt
	;;
    OverlapFirstByLabel)
	run_eval "${BINDIR}/overlap --by-label --first-label 10 parc1.hdr parc2.hdr parc3.hdr > ${tmpdir}/overlap.txt"
	check_result overlap.txt
	;;
    OverlapFirstByLabelNumLabels)
	run_eval "${BINDIR}/overlap --by-label --first-label 10 --num-labels 10 parc1.hdr parc2.hdr parc3.hdr > ${tmpdir}/overlap.txt"
	check_result overlap.txt
	;;
    ProbeXformSimple)
	run_eval "${BINDIR}/probe_xform --probe 180,180,60 --probe 20,20,20 --probe 0,0,0 vol001_mr_t0t1_warp.xform > ${tmpdir}/stdout.txt"
	check_result stdout.txt
	;;
    ProbeXformFwdBwd)
	run_eval "${BINDIR}/probe_xform --probe 180,180,60 --probe 20,20,20 --probe 0,0,0 vol001_mr_t0t1_warp.xform --inverse vol001_mr_t0t1_warp.xform > ${tmpdir}/stdout.txt"
	check_result stdout.txt
	;;
    ProbeXformBwdFwd)
	run_eval "${BINDIR}/probe_xform --probe 180,180,60 --probe 20,20,20 --probe 0,0,0 -- --inverse vol001_mr_t0t1_warp.xform vol001_mr_t0t1_warp.xform > ${tmpdir}/stdout.txt"
	check_result stdout.txt
	;;
    ReformatxNoXform)
	run ${BINDIR}/reformatx --linear -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr
	check_result reformat.img
	;;
    ReformatxLinear)
	run ${BINDIR}/reformatx --linear -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxNearestNeighbor)
	run ${BINDIR}/reformatx --nn -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxPartialVolume)
	run ${BINDIR}/reformatx --pv -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxLinearFwdBwd)
	run ${BINDIR}/reformatx --linear --short -o ${tmpdir}/reformat.hdr --floating vol001_mr_t0.hdr vol001_mr_t0.hdr vol001_mr_t0t1.list --inverse vol001_mr_t0t1.list
	check_result vol001_mr_t0.img
	;;
    ReformatxCubic)
	run ${BINDIR}/reformatx --cubic -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxCubicInverse)
	run ${BINDIR}/reformatx --cubic -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr --inverse vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxSincCosine)
	run ${BINDIR}/reformatx --sinc-cosine -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxSincHamming)
	run ${BINDIR}/reformatx --sinc-hamming -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxSincCosine5)
	run ${BINDIR}/reformatx --sinc-cosine --sinc-window-radius 5 -o ${tmpdir}/reformat.hdr --floating vol001_mr_t1.hdr vol001_mr_t0_crop.hdr vol001_mr_t0t1.list
	check_result reformat.img
	;;
    ReformatxJacobian)
	run ${BINDIR}/reformatx -o ${tmpdir}/jacobian.hdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform --jacobian vol001_mr_t0t1_warp.xform
	check_result jacobian.img
	;;
    ReformatxInverseJacobian)
	run ${BINDIR}/reformatx -o ${tmpdir}/jacobian.hdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform --jacobian --inverse vol001_mr_t0t1_warp.xform
	check_result jacobian.img
	;;
    ReformatxDfieldNrrd)
	run ${BINDIR}/reformatx -o ${tmpdir}/reformat.hdr --floating parc2.hdr parc1.hdr parc1_parc2_dfield.nrrd
	check_result reformat.img
	;;
    RegistrationFromList)
	run ${BINDIR}/registration -v --dofs 0 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    RegistrationWithInitial)
	run ${BINDIR}/registration -v --dofs 0 -o ${tmpdir} --initial vol001_mr_t0t1.list vol001_mr_t0.hdr vol001_mr_t1.hdr
	check_result registration
	;;
    RegistrationWithInitialInverse)
	run ${BINDIR}/registration -v --dofs 0 -o ${tmpdir} --initial vol001_mr_t0t1.list --initial-is-inverse vol001_mr_t0t1.list
	check_result registration
	;;
    RegistrationAutoLevelsRat4)
	run ${BINDIR}/registration -v -i --auto-multi-levels 4 --dofs 6 -o ${tmpdir} rat_fse_erly.hdr rat_fse_late.hdr
	check_result registration
	;;
    RegistrationAutoLevelsRat2)
	run ${BINDIR}/registration -v -i --auto-multi-levels 2 --dofs 6 -o ${tmpdir} rat_fse_erly.hdr rat_fse_late.hdr
	check_result registration
	;;
    RegistrationAutoLevelsRatToRat)
	run ${BINDIR}/registration -v -i --auto-multi-levels 2 --dofs 6,9 --msd --match-histograms -o ${tmpdir} rat_fse_erly.hdr rat2_fse_erly.hdr
	check_result registration
	;;
    RegistrationAutoLevelsRatToRatDeltaFThreshold)
	run ${BINDIR}/registration -v -i --auto-multi-levels 2 --dofs 6,9 --msd --match-histograms --delta-f-threshold 0.01 -o ${tmpdir} rat_fse_erly.hdr rat2_fse_erly.hdr
	check_result registration
	;;
    RegistrationAutoLevelsCt3)
	run ${BINDIR}/registration -q --msd --auto-multi-levels 3 --dofs 6 -o ${tmpdir} pat002_ct.hdr pat002_ct.hdr
	check_result registration
	;;
    ReorientHdrSaToAx)
	run ${BINDIR}/reorient -o RAS phantom_sa.hdr ${tmpdir}/reorient.hdr
	check_result reorient.hdr
	check_result reorient.img
	;;
    ReorientHdrCoToAx)
	run ${BINDIR}/reorient -o RAS phantom_co.hdr ${tmpdir}/reorient.hdr
	check_result reorient.hdr
	check_result reorient.img
	;;
    ReorientHdrAxToSa)
	run ${BINDIR}/reorient -o ASL phantom_ax.hdr ${tmpdir}/reorient.hdr
	check_result reorient.hdr
	check_result reorient.img
	;;
    ReorientHdrCoToSa)
	run ${BINDIR}/reorient -o ASL phantom_co.hdr ${tmpdir}/reorient.hdr
	check_result reorient.hdr
	check_result reorient.img
	;;
    ReorientHdrAxToCo)
	run ${BINDIR}/reorient -o LSA phantom_ax.hdr ${tmpdir}/reorient.hdr
	check_result reorient.hdr
	check_result reorient.img
	;;
    ReorientHdrSaToCo)
	run ${BINDIR}/reorient -o LSA phantom_sa.hdr ${tmpdir}/reorient.hdr
	check_result reorient.hdr
	check_result reorient.img
	;;
    ReorientNrrdToNrrd)
	run ${BINDIR}/reorient vol001_mr_t0_crop.nrrd ${tmpdir}/vol001_mr_t0_crop.nhdr
	check_result vol001_mr_t0_crop.nhdr
	check_result vol001_mr_t0_crop.raw
	;;
    RigidRegistrationMrPet)
	run ${BINDIR}/registration -q -i --dofs 6 -o ${tmpdir} pat001_mr_T1.hdr pat001_pet.hdr
	check_result registration
	;;
    RigidRegistrationMrCt)
	run ${BINDIR}/registration -q -i --dofs 6 -o ${tmpdir} pat002_mr_T2.hdr pat002_ct.hdr
	check_result registration
	;;
    RigidRegistrationCt)
	run ${BINDIR}/registration -q --msd -i --dofs 6 -o ${tmpdir} pat002_ct.hdr pat002_ct.hdr
	check_result registration
	;;
    RigidRegistrationPetMr)
	run ${BINDIR}/registration -q -i --dofs 6 -o ${tmpdir} pat001_pet.hdr pat001_mr_T1.hdr
	check_result registration
	;;
    RigidRegistrationCtMr)
	run ${BINDIR}/registration -q -i --dofs 6 -o ${tmpdir} pat002_ct.hdr pat002_mr_T2.hdr
	check_result registration
	;;
    RigidRegistrationMrPetNoSwap)
	run ${BINDIR}/registration -q --no-switch -i --dofs 6 -o ${tmpdir} pat001_mr_T1.hdr pat001_pet.hdr
	check_result registration
	;;
    RigidRegistrationMrCtNoSwap)
	run ${BINDIR}/registration -q --no-switch -i --dofs 6 -o ${tmpdir} pat002_mr_T2.hdr pat002_ct.hdr
	check_result registration
	;;
    RigidRegistrationPetMrDOF9)
	run ${BINDIR}/registration -q -i --dofs 9 -o ${tmpdir} pat001_pet.hdr pat001_mr_T1.hdr
	check_result registration
	;;
    RigidRegistrationCtMrDOF7)
	run ${BINDIR}/registration -q -i --dofs 7 -o ${tmpdir} pat002_ct.hdr pat002_mr_T2.hdr
	check_result registration
	;;
    RigidRegistrationLabelsDOF69)
	run ${BINDIR}/registration -q --dofs 6 --dofs 9 --class-ref label --class-flt label -o ${tmpdir} parc1.hdr parc2.hdr
	check_result registration
	;;
    RigidRegistrationCrop)
	run ${BINDIR}/registration -v -i -e 2.0 -a 0.125 --sampling 0.25 --crop-index-ref 17,20,0,47,49,12 --crop-index-flt 12,15,0,52,54,12 --dofs 6 -o ${tmpdir} rat_fse_erly.hdr rat_fse_late.hdr
	check_result registration
	;;
    SequenceDefault)
	run_eval "cat numbers.txt | ${BINDIR}/sequence > ${tmpdir}/sequence.txt"
	check_result sequence.txt
	;;
    SequenceFormat)
	run_eval "cat numbers.txt | ${BINDIR}/sequence --format %g > ${tmpdir}/sequence.txt"
	check_result sequence.txt
	;;
    SequenceThresh)
	run_eval "cat numbers.txt | ${BINDIR}/sequence --thresh 1e4 > ${tmpdir}/sequence.txt"	
	check_result sequence.txt
	;;
    SequenceAbs)
	run_eval "cat numbers.txt | ${BINDIR}/sequence --abs > ${tmpdir}/sequence.txt"
	check_result sequence.txt
	;;
    SequenceAbsThresh)
	run_eval "cat numbers.txt | ${BINDIR}/sequence --thresh 1000 --abs > ${tmpdir}/sequence.txt"
	check_result sequence.txt
	;;
    SimilarityGrey)
	run_eval "${BINDIR}/similarity --histogram-file ${tmpdir}/histogram --histogram-text-file ${tmpdir}/histogram.txt rat_fse_erly.hdr rat_fse_late.hdr > ${tmpdir}/similarity.txt"
	check_results histogram histogram.txt similarity.txt
	;;
    SimilarityLabels)
	run_eval "${BINDIR}/similarity --labels --histogram-file ${tmpdir}/histogram --histogram-text-file ${tmpdir}/histogram.txt parc1.hdr parc2.hdr > ${tmpdir}/similarity.txt"
	check_results histogram histogram.txt similarity.txt
	;;
    SimilarityGreyMask)
	run_eval "${BINDIR}/similarity --mask rat_fse_erly.hdr --histogram-file ${tmpdir}/histogram --histogram-text-file ${tmpdir}/histogram.txt rat_fse_erly.hdr rat_fse_late.hdr > ${tmpdir}/similarity.txt"
	check_results histogram histogram.txt similarity.txt
	;;
    SimilarityLabelsMask)
	run_eval "${BINDIR}/similarity --mask parc3_bin.hdr --labels --histogram-file ${tmpdir}/histogram --histogram-text-file ${tmpdir}/histogram.txt parc1.hdr parc2.hdr > ${tmpdir}/similarity.txt"
	check_results histogram histogram.txt similarity.txt
	;;
    SplitAxial)
	run ${BINDIR}/split --output-xform-path ${tmpdir}/split_ax_%1d.xform --axial spgr_3t.hdr ${tmpdir}/split_ax_%1d.hdr
	check_result split_ax_0.img
	check_result split_ax_1.img
	check_result split_ax_0.xform
	check_result split_ax_1.xform
	;;
    SplitAxialNrrd)
	run ${BINDIR}/split --axial --factor 3 spgr_3t.hdr ${tmpdir}/split_ax_%1d.nhdr
	for i in 0 1 2; do
	    check_result split_ax_${i}.nhdr
	    check_result split_ax_${i}.raw
	done
	;;
    SplitSagittal2)
	run ${BINDIR}/split --output-xform-path ${tmpdir}/split_sa_%1d.xform --factor 2 --sagittal spgr_3t.hdr ${tmpdir}/split_sa_%1d.hdr
	check_result split_sa_0.img
	check_result split_sa_1.img
	check_result split_sa_0.xform
	check_result split_sa_1.xform
	;;
    SplitCoronal3)
	run ${BINDIR}/split --output-xform-path ${tmpdir}/split_co_%1d.xform --factor 3 --coronal spgr_3t.hdr ${tmpdir}/split_co_%1d.hdr
	check_result split_co_0.img
	check_result split_co_1.img
	check_result split_co_2.img
	check_result split_co_0.xform
	check_result split_co_1.xform
	check_result split_co_2.xform
	;;
    StatisticsGrey)
	run_eval "${BINDIR}/statistics spgr_3t.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    StatisticsPercentiles)
	run_eval "${BINDIR}/statistics -p 0.1 --percentile 0.5 -p 0.75 spgr_3t.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    StatisticsGreyColumn)
	run_eval "${BINDIR}/statistics -C spgr_3t.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    StatisticsGreyExpNotation)
	run_eval "${BINDIR}/statistics -E spgr_3t.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    StatisticsGreyMask)
	run_eval "${BINDIR}/statistics -m spgr_3t_mask.hdr spgr_3t.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    StatisticsGreyMultiMask)
	run_eval "${BINDIR}/statistics -M spgr_3t_mask.hdr spgr_3t.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    StatisticsLabels)
	run_eval "${BINDIR}/statistics -l parc1.hdr > ${tmpdir}/statistics.txt"
	check_results statistics.txt
	;;
    SymmetryPlane)
	run ${BINDIR}/sympl --sampling 1 --levels 4 --accuracy 0.1 --write-xform ${tmpdir}/xform --sinc --write-subtract ${tmpdir}/subtract.hdr -o ${tmpdir}/parameters cad001_ct.hdr
	check_result parameters
	check_result xform
	check_result subtract.img
	;;
    SymmetryPlaneThresh)
	run ${BINDIR}/sympl --sampling 1 --levels 4 --accuracy 0.1 --min-value -224 --max-value 176 --cubic --write-aligned ${tmpdir}/aligned.hdr -o ${tmpdir}/parameters cad001_ct.hdr
	check_result parameters
	check_result aligned.img
	;;
    TTestDefault)
	run ${BINDIR}/ttest -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestOneSided)
	run ${BINDIR}/ttest -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestSymmetric)
	run ${BINDIR}/ttest -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr --symmetric jacobian-01.nii jacobian-02.nii jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestLog)
	run ${BINDIR}/ttest --log -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestAbsLog)
	run ${BINDIR}/ttest --abs --log -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestInvert)
	run ${BINDIR}/ttest --invert -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestPaired)
	run ${BINDIR}/ttest --paired -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    TTestCrossCorrelation)
	run ${BINDIR}/ttest --cross-correlation -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/pvals.hdr jacobian-01.nii jacobian-02.nii jacobian-03.nii -- jacobian-04.nii jacobian-03.nii jacobian-02.nii 
	check_results ttest.img pvals.img
	;;
    TTestZScores)
	run ${BINDIR}/ttest --zscores -o ${tmpdir}/zscores.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results zscores.img
	;;
    TTestMask)
	run ${BINDIR}/ttest --mask jacobian-mask.nii -o ${tmpdir}/ttest.hdr --tstats-file ${tmpdir}/tstats.hdr jacobian-01.nii jacobian-02.nii -- jacobian-03.nii jacobian-04.nii
	check_results ttest.img tstats.img
	;;
    UnsplitHdrAx)
	run ${BINDIR}/unsplit --axial -o ${tmpdir}/unsplit.hdr split_ax_0.hdr split_ax_1.hdr
	check_result unsplit.hdr
	check_result unsplit.img
	;;
    UnsplitHdrSa)
	run ${BINDIR}/unsplit --sagittal -o ${tmpdir}/unsplit.hdr split_sa_0.hdr split_sa_1.hdr
	check_result unsplit.hdr
	check_result unsplit.img
	;;
    UnsplitHdrCo)
	run ${BINDIR}/unsplit --coronal -o ${tmpdir}/unsplit.hdr split_co_0.hdr split_co_1.hdr split_co_2.hdr
	check_result unsplit.hdr
	check_result unsplit.img
	;;
    UnsplitHdrNrrdAx)
	run ${BINDIR}/unsplit --axial -o ${tmpdir}/unsplit.nhdr split_ax_0.hdr split_ax_1.hdr
	check_result unsplit.nhdr
	check_result unsplit.raw
	;;
    UnsplitHdrNrrdSa)
	run ${BINDIR}/unsplit --sagittal -o ${tmpdir}/unsplit.nhdr split_sa_0.hdr split_sa_1.hdr
	check_result unsplit.nhdr
	check_result unsplit.raw
	;;
    UnsplitHdrNrrdCo)
	run ${BINDIR}/unsplit --coronal -o ${tmpdir}/unsplit.nhdr split_co_0.hdr split_co_1.hdr split_co_2.hdr
	check_result unsplit.nhdr
	check_result unsplit.raw
	;;
    UnsplitNrrdNrrd)
	run ${BINDIR}/unsplit --axial -o ${tmpdir}/unsplit.nhdr split_ax_0.nhdr split_ax_1.nhdr split_ax_2.nhdr
	check_result unsplit.nhdr
	check_result unsplit.raw
	;;
    VolumeInjection)
	run ${BINDIR}/volume_injection --recon-grid-path spgr_3t.hdr -o ${tmpdir}/injection.hdr --gauss-sigma 0.5 --radius 2 split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	;;
    VolumeInjectionIsotropic)
	run ${BINDIR}/volume_injection --recon-grid-path spgr_3t.hdr -o ${tmpdir}/injection.hdr --gauss-sigma 0.5 --radius 2 --isotropic-injection split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	;;
    VolumeInjectionNoXform)
	run ${BINDIR}/volume_injection -o ${tmpdir}/injection.hdr --gauss-sigma 0.5 --radius 2 --exclude-first-image spgr_3t.hdr -- split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	;;
    VolumeInjectionNoXformIsotropic)
	run ${BINDIR}/volume_injection -o ${tmpdir}/injection.hdr --gauss-sigma 0.5 --radius 2 --isotropic-injection --exclude-first-image spgr_3t.hdr -- split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	;;
    VolumeReconstructionFourthOrder)
	run ${BINDIR}/volume_reconstruction --recon-grid-path spgr_3t.hdr -o ${tmpdir}/reconstruction.hdr --fourth-order-error --num-iterations 2 --gauss-sigma 0.5 --radius 2 --isotropic-injection --write-splatted-image ${tmpdir}/injection.hdr split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	check_result reconstruction.img
	;;
    VolumeReconstructionCubic)
	run ${BINDIR}/volume_reconstruction --recon-grid-path spgr_3t.hdr -o ${tmpdir}/reconstruction.hdr --cubic --num-iterations 2 --gauss-sigma 0.5 --radius 2 --isotropic-injection --write-splatted-image ${tmpdir}/injection.hdr split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	check_result reconstruction.img
	;;
    VolumeReconstructionNoXform)
	run ${BINDIR}/volume_reconstruction -o ${tmpdir}/reconstruction.hdr --linear --num-iterations 2 --gauss-sigma 0.5 --radius 2 --isotropic-injection --write-splatted-image ${tmpdir}/injection.hdr --exclude-first-image spgr_3t.hdr -- split_ax_0.hdr split_ax_01.xform split_ax_1.hdr
	check_result injection.img
	check_result reconstruction.img
	;;
    WarpSingleLevel)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 160 --accuracy 1 --no-adaptive-fix -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpSingleLevelExact)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 180 --exact-spacing --accuracy 1 --sampling 3 --no-adaptive-fix -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpInverseConsistentCC)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 80 --accuracy 1 --ncc --ic-weight 1e-2 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpMultiLevel)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 160 --accuracy 1 --refine 1 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpMultiLevelMatchHistograms)
	run ${BINDIR}/warp -q --exploration 8 --match-histograms --msd --grid-spacing 160 --accuracy 1 --refine 1 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpMultiLevelDeltaFThreshold)
	run ${BINDIR}/warp -q --exploration 8 --delta-f-threshold 0.01 --msd --grid-spacing 160 --accuracy 1 --refine 1 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpMultiLevelExact)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 160 --exact-spacing --accuracy 1 --refine 1 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpDelayRefine)
	run ${BINDIR}/warp -q --exploration 12 --grid-spacing 160 --accuracy 2 --refine 1 --delay-refine --sampling 6 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpEnergy)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 160 --accuracy 1 --refine 1 --energy-weight 1e-1 -o ${tmpdir} vol001_mr_t0t1.list
	check_result registration
	;;
    WarpJacobian)
	export CMTK_NUM_THREADS=1
	run ${BINDIR}/warp -q --exploration 12 --grid-spacing 160 --accuracy 2 --refine 1 --jacobian-weight 1e-1 --sampling 12 --omit-original-data -o ${tmpdir} vol001_mr_t0t1.list
	unset CMTK_NUM_THREADS
	check_result registration
	;;
    WarpLabels)
	run ${BINDIR}/warp -q --exploration 8 --grid-spacing 90 --accuracy 1 --refine 1 --class-ref label --class-flt label -o ${tmpdir} --initial parc1_parc2_9dof.xform parc1.hdr parc2.hdr
	check_result registration
	;;
    Xform2dfieldWarpNrrd)
	run ${BINDIR}/xform2dfield -v ${tmpdir}/dfield.nhdr vol001_mr_t0_crop.hdr vol001_mr_t0t1_warp.xform
	check_result dfield.nhdr
	check_result dfield.raw
	;;
    Xform2dfieldAffineNrrd)
	run ${BINDIR}/xform2dfield -v ${tmpdir}/dfield.nhdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform
	check_result dfield.nhdr
	check_result dfield.raw
	;;
    Xform2dfieldDownsampleXYZNrrd)
	run ${BINDIR}/xform2dfield -v --downsample 4,4,2 ${tmpdir}/dfield.nhdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform
	check_result dfield.nhdr
	check_result dfield.raw
	;;
    Xform2dfieldDownsampleXNrrd)
	run ${BINDIR}/xform2dfield -v --downsample 4 ${tmpdir}/dfield.nhdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform
	check_result dfield.nhdr
	check_result dfield.raw
	;;
    Xform2dfieldConcatNrrd)
	run ${BINDIR}/xform2dfield -v ${tmpdir}/dfield.nhdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform vol001_mr_t0t1_warp.xform
	check_result dfield.nhdr
	check_result dfield.raw
	;;
    Xform2dfieldInverseNrrd)
	run ${BINDIR}/xform2dfield -v ${tmpdir}/dfield.nhdr vol001_mr_t0_crop.hdr vol001_mr_t0_crop.xform --inverse vol001_mr_t0t1_warp.xform
	check_result dfield.nhdr
	check_result dfield.raw
	;;
    xml_film)
	run_eval "${BINDIR}/film --xml | sed '/<version>/{ N; s/^.*$/<version>/ }' > ${tmpdir}/film.xml"
	check_result film.xml
	;;
    xml_levelset)
        run_eval "${BINDIR}/levelset --xml | sed '/<version>/{ N; s/^.*$/<version>/ }' > ${tmpdir}/levelset.xml"
	check_result levelset.xml
	;;
    xml_mrbias)
	run_eval "${BINDIR}/mrbias --xml | sed '/<version>/{ N; s/^.*$/<version>/ }' > ${tmpdir}/mrbias.xml"
	check_result mrbias.xml
	;;
    xml_registration)
	run_eval "${BINDIR}/registration --xml | sed '/<version>/{ N; s/^.*$/<version>/ }' > ${tmpdir}/registration.xml"
	check_result registration.xml
	;;
    wiki_film)
	run_eval "${BINDIR}/film --wiki > ${tmpdir}/film.wiki"
	check_result film.wiki
	;;
    wiki_levelset)
	run_eval "${BINDIR}/levelset --wiki > ${tmpdir}/levelset.wiki"
	check_result levelset.wiki
	;;
    wiki_mrbias)
	run_eval "${BINDIR}/mrbias --wiki > ${tmpdir}/mrbias.wiki"
	check_result mrbias.wiki
	;;
    wiki_registration)
	run_eval "${BINDIR}/registration --wiki > ${tmpdir}/registration.wiki"
	check_result registration.wiki
	;;
    *)
	exit 2
	;;
esac

if [ "${tmpdir}" != "" ]; then
    rm -rf ${tmpdir}
fi
