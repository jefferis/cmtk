##
##  Copyright 2004-2012 SRI International
##
##  Copyright 1997-2011 Torsten Rohlfing
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

FIND_PROGRAM(SH_PATH sh PATHS C:/cygwin/bin)
SET(testDriver ${SH_PATH} ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh)
SET(testList "")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.sh.in ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh @ONLY)
# On Windows, remove "\r" (CR) from line ends, which would make Cygwin's sh choke
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  EXECUTE_PROCESS(COMMAND c:\\cygwin\\bin\\sed -i.bak s/\\r//g ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# ==========================================
# Tests for "avg_adm" tool
LIST(APPEND testList
  avg_admDefault2
  avg_admDefault3
  avg_admJacobianFloat
  avg_admLabels
  avg_admPaddOutZeroNN
  avg_admNoReferenceModelCubic
  avg_admNoScaleModelAutoScale
  avg_admWithAffineNoRefData)

# ==========================================
# Tests for "dbtool"
IF(CMTK_USE_SQLITE)
  LIST(APPEND testList
    dbtool_AddImages
    dbtool_AddImages2
    dbtool_ListSpace
    dbtool_GetXform1
    dbtool_GetXform2
    dbtool_GetXform3
    dbtool_GetXform4
    dbtool_GetXform5
    dbtool_GetXform6)
ENDIF(CMTK_USE_SQLITE)  

# ==========================================
# Tests for "describe" tool
LIST(APPEND testList
  describeMountPoints
  describeMountPointsMulti
  describeMountPointsPrefix
  describeMountPointsPrefixInvalid
  describeBZip2a
  describeBZip2b
  describeLZMAa 
  describeLZMAb
  describeXZa 
  describeXZb
  describeEmpty
  describeDICOM
  describeMosaicDICOM
  describeDICOMZ
  describeVanderbilt
  describeMR1 
  describeMR2
  describeMR3
  describeMR4
  describeHuman
  describeMRBiorad
  describeMRBioradGz
  describeNiftiDetached348
  describeEmbedAnalyze
  describeEmbedNifti)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    describeMRNrrd1 
    describeMRNrrd2
    describeEmbedNrrd)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "destripe" tool
LIST(APPEND testList
  destripeDefault
  destripeKernelFloat)

# ==========================================
# Tests for "detect_adni_phantom" tool
IF(CMTK_USE_FFTW)
  LIST(APPEND testList 
    detect_adni_phantomDefault)
ENDIF(CMTK_USE_FFTW)

# ==========================================
# Tests for "detect_spheres_matched_filter" tool
IF(CMTK_USE_FFTW)
  LIST(APPEND testList 
    detect_spheres_matched_filter
    detect_spheres_matched_filterNormalized)
ENDIF(CMTK_USE_FFTW)

# ==========================================
# Tests for "epiunwarp" tool
IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList epiunwarp epiunwarpInitShift)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "fib2image" tool
LIST(APPEND testList fib2image)

# ==========================================
# Tests for "fibxform" tool
LIST(APPEND testList fibxform fibxformSourceTarget)

# ==========================================
# Tests for "groupwise_init" tool
LIST(APPEND testList
  groupwise_initCenterFOV
  groupwise_initCentersOfMass   
  groupwise_initCentersOfMassScale
  groupwise_initCentersOfMassTemplate)

# ==========================================
# Tests for "groupwise_affine" tool
LIST(APPEND testList
  groupwise_affineFromInit
  groupwise_affineMatchHistograms
  groupwise_affineFromInitSampling
  groupwise_affineBackground 
  groupwise_affineUseTemplate
  groupwise_affineZeroSumSmooth
  groupwise_affineRMIFromInit
  groupwise_affineRMIFromInitDeltaF
  groupwise_affineRMIFromInitSampling 
  groupwise_affineRMIBackground 
  groupwise_affineRMIZeroSumSmooth)

# ==========================================
# Tests for "groupwise_warp" tool
LIST(APPEND testList
  groupwise_warpFromInit
  groupwise_warpFromInitProtect
  groupwise_warpFitFromInit
  groupwise_warpFromInitZeroSum
  groupwise_warpUseTemplate
  groupwise_warpUseTemplateMatchHistograms
  groupwise_warpMatchHistograms
  groupwise_warpFromInitNoBG
  groupwise_warpUseTemplateNoBG
  groupwise_warpRMIFromInit
  groupwise_warpRMIFitFromInit
  groupwise_warpRMIFromInitZeroSum
  )

# ==========================================
# Tests for xml, wiki, man, and help markups
IF(CMTK_USE_SQLITE)
  FOREACH(t film levelset mrbias registration)
    LIST(APPEND testList xml_${t} wiki_${t} man_${t})
  ENDFOREACH(t film levelset mrbias registration)
  LIST(APPEND xml_gmm)
ENDIF(CMTK_USE_SQLITE)

FOREACH(t film)
  LIST(APPEND testList help_${t} help_all_${t})
ENDFOREACH(t film)

# ==========================================
# Tests for "concat_affine" tool
LIST(APPEND testList
  concat_affineABA
  concat_affineABAInvert
  concat_affineAB1A	     
  concat_affineAA1
  concat_affineA1A)

# ==========================================
# Tests for "average_images" tool
LIST(APPEND testList
  average_imagesMean 
  average_imagesMeanNormPadd 
  average_imagesMeanAbsLog 
  average_imagesVariance
  average_imagesStDev
  average_imagesEntropy
  average_imagesZScore)

# ==========================================
# Tests for "convertx" tool
LIST(APPEND testList
  convertxType
  convertxBinarize
  convertxBinarizeOtsu
  convertxBinarizeOtsuNBins
  convertxPruneHighLow
  convertxPruneHigh
  convertxPruneLow
  convertxBoundaryMap
  convertxBoundaryMapMultiValue
  convertxConnectedComponents
  convertxCropThresholdWriteRegion
  convertxCropRegion1
  convertxCropRegion2
  convertxDownsample
  convertxDownsampleSelect
  convertxDistanceUnsigned
  convertxDistanceSigned
  convertxDistanceSigned2
  convertxFlipX
  convertxFlipYZ
  convertxErodeDilateErode
  convertxDilateErodeDilate
  convertxGaussianFilterSigma
  convertxGaussianFilterFWHM
  convertxHistogramEqualization
  convertxHistogramEqualizationNBins
  convertxLaplace
  convertxMapValues
  convertxMapValues2
  convertxMapValues3
  convertxMapValuesOnly
  convertxMapValuesOnly2
  convertxMask
  convertxMaskInverse
  convertxMedianFilter1
  convertxMedianFilter2
  convertxMedianFilterXYZ
  convertxMeanFilter
  convertxFastMean0Filter
  convertxFastMean1Filter
  convertxFastMean2Filter
  convertxFastMean2PadFilter
  convertxVarianceFilter
  convertxFastVariance0Filter
  convertxFastVariance1Filter
  convertxFastVariance2Filter
  convertxFastVariance2PadFilter
  convertxThirdMomentFilter
  convertxStandardDeviationFilter
  convertxSmoothnessFilter
  convertxMedialSkeletonFilter1D
  convertxMedialSkeletonFilter2D
  convertxNiftiToAnalyze
  convertxNiftiToMetaImage
  convertxAnalyzeToNifti
  convertxAnalyzeToNiftiRAS
  convertxNiftiDetachedToNifti
  convertxAnalyzeToNiftiDetached 
  convertxUncompressedNIFTI1
  convertxUncompressedNIFTI2
  convertxUncompressedNIFTI3
  convertxUncompressedNIFTI4
  convertxThresholdBelow
  convertxThresholdAbove
  convertxThresholdBelowToPadding
  convertxThresholdAboveToPadding
  convertxScaleToRange
  convertxRevert
  convertxReplacePadding    
  convertxReplaceInfNaN)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList convertxDownsampleNrrd convertxDownsampleSelectNrrd)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "dcm2image" tool
IF(CMTK_USE_DCMTK)
  LIST(APPEND testList
    dcm2image 
    dcm2imageSubs
    dcm2imageZ
    dcm2imageEmbedPatientNameAnalyze
    dcm2imageEmbedPatientNameNifti
    dcm2imageEmbedSeriesDescriptionNifti
    dcm2imageEmbedStudyIDDateNifti
    dcm2imageEmbedDefaultNifti
    dcm2imageMosaic
    dcm2imageDiffusionGEXML
    dcm2imageExclude
    dcm2imageInclude)
  
  IF(CMTK_BUILD_NRRD)
    LIST(APPEND testList
      dcm2imageNrrd
      dcm2imageEmbedPatientNameNrrd)
  ENDIF(CMTK_BUILD_NRRD)
ENDIF(CMTK_USE_DCMTK)

# ==========================================
# Tests for "dof2mat" tool
LIST(APPEND testList dof2mat dof2matTranspose)

# ==========================================
# Tests for "filter" tool
LIST(APPEND testList
  filterGaussian
  filterGaussianSmallKernel
  filterGaussianNoFilter
  filterRohlfing)

# ==========================================
# Tests for "film" tool
LIST(APPEND testList
  filmCubic
  filmFourthOrder
  filmMSDLinearNoTrunc
  filmMISincRefSPGR)


# ==========================================
# Tests for "fit_affine_xform_landmarks" tool
LIST(APPEND testList fit_affine_xform_landmarksRigid fit_affine_xform_landmarksScales fit_affine_xform_landmarksShear
  fit_affine_xform_landmarksRigidRigid fit_affine_xform_landmarksRigidRigidFlip fit_affine_xform_landmarksRigidScales fit_affine_xform_landmarksRigidShear)

# ==========================================
# Tests for "fit_affine_xform" tool
LIST(APPEND testList fit_affine_xformFFD)
IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList fit_affine_xformDField)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "fit_spline_xform" tool
IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    fit_spline_xformSpacing
    fit_spline_xformSpacingMultiLevel
    fit_spline_xformDims
    fit_spline_xformDimsRelative
    fit_spline_xformDimsMultiLevel
    fit_spline_xformDimsMultiLevelSafe
    fit_spline_xformDimsWithAffine)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "glm" tool
LIST(APPEND testList
  glmDefault glmNormalize glmExp glmNoConstant glmIgnore glmSelect glmCrop)

# ==========================================
# Tests for "gmm" tool
LIST(APPEND testList gmmDefault gmmAlternative)

# ==========================================
# Tests for "gregxform" tool
LIST(APPEND testList
  gregxformFordwardBackward 
  gregxformAffine 
  gregxformAffineFromWarp
  gregxformAffineFromWarpFwdBwd)

# ==========================================
# Tests for "hausdorff" tool
LIST(APPEND testList
  hausdorffBinary12
  hausdorffBinary21
  hausdorffBinaryForced)

# ==========================================
# Tests for "histogram" tool
LIST(APPEND testList
  histogram 
  histogramNorm
  histogramBinsMinMax
  histogramBinsMinMaxTrunc
  histogramMask)

# ==========================================
# Tests for "imagemath" tool
LIST(APPEND testList
  imagemathAbsPop
  imagemathSqrSqrt
  imagemathLogExp
  imagemathDupAddMulDiv
  imagemathThreshAboveBelow
  imagemathAverage 
  imagemathVariance
  imagemathMaskAverage 
  imagemathXor
  imagemathDupFillMaxValue
  imagemathAnd
  imagemathContractLabels
  imagemathMaxIndex
  imagemathMaxValue
  imagemathStackEntropyLabels
  imagemathVote
  imagemathProduct
  imagemathSTAPLE
  imagemathMultiClassSTAPLE 
  imagemathMultiClassDisputedSTAPLE 
  imagemathCombinePCA
  imagemathT2
  imagemathAtan2
  imagemathLogOddsAdd
  imagemathLogOddsAdd2
  imagemathMatchMeanSDev
  imagemathMatchMeanSDev3
  imagemathMatchHistograms 
  imagemathMatchHistogramsPadding
  imagemathMatchHistogramsPadding2
  imagemathMatchHistogramsPaddingUnset)

# ==========================================
# Tests for "jidb" tool
LIST(APPEND testList
  jidbGaussian
  jidbBoxFourthOrder
  jidbGaussianScale
  jidbMIRefSPGR)

# ==========================================
# Tests for "levelset" tool
LIST(APPEND testList
  levelsetDefault
  levelsetScaleInitial
  levelsetBinarizeFastWideBinary
  )

# ==========================================
# Tests for "lmsba" tool
LIST(APPEND testList lmsbaRadius3 lmsbaRadius3Outliers lmsbaRadius3OutliersGlobal lmsbaRadius3Search1)

# ==========================================
# Tests for "lsba" tool
LIST(APPEND testList lsbaRadius3 lsbaRadius3Outliers lsbaRadius3OutliersGlobal lsbaRadius3Search1)

# ==========================================
# Tests for "lvote" tool
LIST(APPEND testList lvoteDefault lvoteGlobalWeights lvoteRadius3 lvoteOutliersGlobal)

# ==========================================
# Tests for "make_initial_affine" tool
LIST(APPEND testList
  make_initial_affineCenterOfMass
  make_initial_affinePrincipalAxes1
  make_initial_affinePrincipalAxes2
  make_initial_affinePrincipalAxes3
  make_initial_affinePrincipalAxes4
  make_initial_affinePrincipalAxes5
  make_initial_affinePrincipalAxes6)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    make_initial_affineDirectionVectorsNrrdAxSa
    make_initial_affineDirectionVectorsNrrdAxCo
    make_initial_affineDirectionVectorsNrrdSaCo
    make_initial_affineDirectionVectorsNrrdAxSaNative
    make_initial_affineDirectionVectorsNrrdAxCoNative
    make_initial_affineDirectionVectorsNrrdSaCoNative)
ENDIF(CMTK_BUILD_NRRD)

IF(CMTK_USE_SQLITE)
  LIST(APPEND testList
    make_initial_affineDB)
ENDIF(CMTK_USE_SQLITE)  

# ==========================================
# Tests for "mat2dof" tool
LIST(APPEND testList
  mat2dof1
  mat2dof2
  mat2dofFile)

# ==========================================
# Multi-channel registration requires hash map to build
IF(HAVE_STL_HASH_MAP)
  # ==========================================
  # Tests for "mcaffine" tool
  LIST(APPEND testList
    mcaffine1 mcaffine2 mcaffine3 mcaffine4)
  
  # ==========================================
  # Tests for "mcwarp" tool
  LIST(APPEND testList
    mcwarp1 mcwarp2 mcwarp3 mcwarp4 mcwarp5 mcwarp6)

ENDIF(HAVE_STL_HASH_MAP)

# ==========================================
# Tests for "mip" tool
LIST(APPEND testList
  mipDefault
  mipBlackWhite
  mipAdjust
  mipSum
  mip16Bit)

# ==========================================
# Tests for "mk_adni_phantom" tool
LIST(APPEND testList mk_adni_phantom2mm mk_adni_phantom2mmLabels)

# ==========================================
# Tests for "mk_analyze_hdr" tool
LIST(APPEND testList
  mk_analyze_hdrDefault
  mk_analyze_hdrImport
  mk_analyze_hdrLittleEndian
  mk_analyze_hdrBigEndian)

# ==========================================
# Tests for "mk_nifti_hdr" tool
LIST(APPEND testList
  mk_nifti_hdrDefault
  mk_nifti_hdrDescription
  mk_nifti_hdrImport
  mk_nifti_hdrDefaultAttached
  mk_nifti_hdrDescriptionAttached
  mk_nifti_hdrImportAttached)

# ==========================================
# Tests for "mk_phantom_3d" tool
LIST(APPEND testList
  mk_phantom_3dBoxIndexed 
  mk_phantom_3dBoxIndexedRange
  mk_phantom_3dBoxAbsolute
  mk_phantom_3dBoxRelative
  mk_phantom_3dSphereIndexed
  mk_phantom_3dSphereAbsolute
  mk_phantom_3dSphereAbsolute2
  mk_phantom_3dSphereRelative
  mk_phantom_3dSphereRelative2
  mk_phantom_3dBoxSphere 
  mk_phantom_3dImport
  mk_phantom_3dImportGrid)

# ==========================================
# Tests for "reorient" tool
LIST(APPEND testList
  reorientHdrSaToAx
  reorientHdrCoToAx
  reorientHdrAxToSa
  reorientHdrCoToSa
  reorientHdrAxToCo
  reorientHdrSaToCo)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    reorientNrrdToNrrd
    reorientNrrdToNrrdRAS
    reorientNrrdToNrrdSpaceLPS)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "split" tool
LIST(APPEND testList splitAxial splitSagittal2 splitCoronal3)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList splitAxialNrrd)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "overlap" tool
LIST(APPEND testList
  overlap
  overlapNumLabels
  overlapByLabel
  overlapFirst
  overlapFirstByLabel
  overlapFirstByLabelNumLabels)

# ==========================================
# Tests for "probe" tool
LIST(APPEND testList
  probeIndexed
  probeAbsolute
  probeRelative
  probePhysical
  probeIndexedLinear
  probeIndexedCubic
  probeIndexedSinc)

# ==========================================
# Tests for "registration" tool
LIST(APPEND testList
  registrationAffineMrMrMSD
  registrationRigidMrPet 
  registrationRigidMrCt
  registrationRigidCt
  registrationRigidPetMr
  registrationRigidCtMr
  registrationRigidMrPetNoSwap
  registrationRigidMrCtNoSwap
  registrationRigidPetMrDOF9
  registrationRigidCtMrDOF7
  registrationRigidLabelsDOF69
  registrationRigidCrop
  registrationAutoLevelsRat4
  registrationAutoLevelsRat2
  registrationAutoLevelsRatToRat
  registrationAutoLevelsRatToRatDeltaFThreshold
  registrationAutoLevelsCt3
  registrationFromList
  registrationWithInitial
  registrationWithInitialInverse)

# ==========================================
# Tests for "registrationx" tool
LIST(APPEND testList
  registrationxAffineMrMrMSD
  registrationxShearNoScaleMrMrMSD
  registrationxFromList
  registrationxWithInitial
  registrationxWithInitialInverse
  registrationxAutoLevelsRat4
  registrationxAutoLevelsRat4XY
  registrationxAutoLevelsRat4YZ
  registrationxAutoLevelsRat4XZ
  registrationxAutoLevelsRat4Symmetric
  registrationxAutoLevelsRat4NONE
  registrationxAutoLevelsRat4FOV
  registrationxAutoLevelsRat4COM
  registrationxAutoLevelsRat4PAX
  registrationxAutoLevelsRat2
  registrationxAutoLevelsRat2Cubic
  registrationxAutoLevelsRat2Sinc
  registrationxAutoLevelsRat2NN
  registrationxAutoLevelsRatToRat
  registrationxAutoLevelsRatToRatDeltaFThreshold
  registrationxAutoLevelsRatToRatRMS
  registrationxAutoLevelsRatToRatNCC
  registrationxAutoLevelsRatToRatCR
  registrationxAutoLevelsRatToRatMI
  registrationxAutoLevelsLabelsNN
  registrationxAutoLevelsLabelsPV
  registrationxAutoLevelsCt3)

IF(CMTK_USE_SQLITE)
  LIST(APPEND testList
    registrationxFromListDB
    registrationxWithInitialDB
    registrationxWithInitialInverseDB
    )
ENDIF(CMTK_USE_SQLITE)  

# ==========================================
# Tests for "mrbias" tool
LIST(APPEND testList
  mrbiasMulIncremental
  mrbiasMulAutoThresh
  mrbiasMulLogIntensity
  mrbiasAddMulMask)

IF(CMTK_USE_CUDA)
  LIST(APPEND testList
    mrbiasMulIncrementalCUDA
    mrbiasMulLogIntensityCUDA
    mrbiasAddMulMaskCUDA)
ENDIF(CMTK_USE_CUDA)

# ==========================================
# Tests for "reformatx" tool
LIST(APPEND testList
  reformatxNoXform
  reformatxLinear 
  ##		   reformatxLinearFwdBwd
  reformatxNearestNeighbor
  reformatxPartialVolume
  reformatxCubic
  reformatxCubicInverse
  reformatxSincCosine
  reformatxSincCosine5
  reformatxMassPreserving
  reformatxJacobian
  reformatxInverseJacobian
  reformatxSincHamming
  reformatxTargetGrid
  reformatxTargetGridAnalyze
  reformatxTargetGridOffset)
IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    reformatxDfieldNrrd
    reformatxDfieldNrrdJacobian)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "sba" tool
LIST(APPEND testList sbaDefault sbaOutliers sbaOutliers2)

# ==========================================
# Tests for "sbai" tool
LIST(APPEND testList sbaiDefault)

# ==========================================
# Tests for "sequence" tool
LIST(APPEND testList sequenceDefault sequenceFormat sequenceThresh sequenceAbs sequenceAbsThresh)

# ==========================================
# Tests for "similarity" tool
LIST(APPEND testList similarityGrey similarityLabels similarityGreyMask similarityLabelsMask)

# ==========================================
# Tests for "streamxform" tool
LIST(APPEND testList streamxformFordwardBackward streamxformFordwardBackwardTolerance streamxformAffine streamxformAffineFromTo)

# ==========================================
# Tests for "sympl" tool
LIST(APPEND testList symmetry_plane symmetry_planeThresh)

# ==========================================
# Tests for "symplx" tool
LIST(APPEND testList
  symplx_Default
  symplx_Thresh
  symplx_FixOffset)

IF(CMTK_USE_CUDA)
  LIST(APPEND testList
    symplx_DefaultCUDA
    symplx_ThreshCUDA
    symplx_FixOffsetCUDA)
ENDIF(CMTK_USE_CUDA)

# ==========================================
# Tests for "statistics" tool
LIST(APPEND testList
  statisticsGrey 
  statisticsPercentiles
  statisticsGreyColumn 
  statisticsGreyExpNotation
  statisticsGreyMask 
  statisticsGreyMultiMask
  statisticsMaskMismatch
  statisticsLabels
  statisticsLabelsAllUpToHi
  statisticsLabelsAllUpToLo)

# ==========================================
# Tests for "ttest" tool
LIST(APPEND testList ttestDefault ttestLog ttestAbsLog ttestInvert
  ttestPaired ttestCrossCorrelation ttestZScores
  ttestMask ttestOneSided)

# ==========================================
# Tests for "unsplit" tool
LIST(APPEND testList unsplitHdrAx unsplitHdrSa unsplitHdrCo)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList unsplitHdrNrrdAx unsplitHdrNrrdSa unsplitHdrNrrdCo unsplitNrrdNrrd)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "unwarp_image_phantom" tool
LIST(APPEND testList
  unwarp_image_phantomSpacing
  unwarp_image_phantomSpacingLevels
  unwarp_image_phantomSpacingDirect
  unwarp_image_phantomDims
  unwarp_image_phantomDimsLevels)

# ==========================================
# Tests for "volume_injection" tool
LIST(APPEND testList
  volume_injection 
  volume_injectionReconGrid
  volume_injectionReconGridOffset
  volume_injectionIsotropic 
  volume_injectionNoXform 
  volume_injectionNoXformIsotropic)

# ==========================================
# Tests for "volume_reconstruction" tool
LIST(APPEND testList
  volume_reconstructionFourthOrder 
  volume_reconstructionCubic
  volume_reconstructionNoXform
  volume_reconstructionBoxPSF
  volume_reconstructionGaussianPSF)

# ==========================================
# Tests for "warp" tool
LIST(APPEND testList
  warpSingleLevel
  warpSingleLevelExact
  warpInverseConsistentCC
  warpMultiLevel
  warpMultiLevelMatchHistograms
  warpMultiLevelDeltaFThreshold
  warpMultiLevelExact
  warpDelayRefine
  warpEnergy
  warpEnergyRelax
  warpRigidity
  warpJacobian
  warpLabels)

# ==========================================
# Tests for "warpx" tool
LIST(APPEND testList
  warpxSingleLevel
  warpxSingleLevelExact
  warpxInverseConsistentCC
  warpxMultiLevel
  warpxMultiLevelMatchHistograms
  warpxMultiLevelDeltaFThreshold
  warpxMultiLevelExact
  warpxDelayRefine
  warpxEnergy
  warpxEnergyRelax
  warpxJacobian
  warpxJacobianUnfold
  warpxLabels)

# ==========================================
# Tests for "xform2dfield" tool
IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    xform2dfieldWarpNrrd
    xform2dfieldAffineNrrd
    xform2dfieldDownsampleXYZNrrd
    xform2dfieldDownsampleXNrrd
    xform2dfieldConcatNrrd
    xform2dfieldInverseNrrd)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "xform2scalar" tool
LIST(APPEND testList
  xform2scalarAffine 
  xform2scalarAffineDoubleY 
  xform2scalarWarp 
  xform2scalarWarpInverseError
  xform2scalarWarpOnly)

# ==========================================
# Tests for "vtkxform" tool
LIST(APPEND testList
  vtkxform
  vtkxformInverse)

IF(CMTK_BUILD_NRRD)
  LIST(APPEND testList
    xform2scalarDfield)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Set up all tests
FOREACH(testName ${testList})
  IF(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(NAME ${testName} COMMAND ${testDriver} ${testName} $<CONFIGURATION> ${MEMORYCHECK_COMMAND})
  ELSE(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(NAME ${testName} COMMAND ${testDriver} ${testName} $<CONFIGURATION>)
  ENDIF(CMTK_TESTING_MEMORYCHECK)
ENDFOREACH(testName ${testList})

SET_TESTS_PROPERTIES(groupwise_warpFromInit PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(groupwise_warpFitFromInit PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(groupwise_warpFromInitZeroSum PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(groupwise_warpUseTemplate PROPERTIES TIMEOUT 6400)
SET_TESTS_PROPERTIES(groupwise_warpFromInitNoBG PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(groupwise_warpUseTemplateNoBG PROPERTIES TIMEOUT 6400)

SET_TESTS_PROPERTIES(groupwise_warpRMIFromInit groupwise_warpRMIFitFromInit groupwise_warpRMIFromInitZeroSum PROPERTIES TIMEOUT 4800)

# ==========================================
# Set up test properties, now that the tests 
# have been defined.

SET_TESTS_PROPERTIES(avg_admLabels filmFourthOrder jidbGaussian jidbGaussianScale volume_reconstructionGaussianPSF volume_reconstructionNoXform PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(filmCubic PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(filmMSDLinearNoTrunc PROPERTIES TIMEOUT 3600)
SET_TESTS_PROPERTIES(filmMISincRefSPGR PROPERTIES TIMEOUT 3600)

SET_TESTS_PROPERTIES(jidbMIRefSPGR PROPERTIES TIMEOUT 3600)

SET_TESTS_PROPERTIES(registrationxShearNoScaleMrMrMSD PROPERTIES TIMEOUT 900)

IF(HAVE_STL_HASH_MAP)
  SET_TESTS_PROPERTIES(mcaffine2 PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(mcwarp1 mcwarp2 mcwarp3 mcwarp4 mcwarp5 mcwarp6 PROPERTIES TIMEOUT 1800)
ENDIF(HAVE_STL_HASH_MAP)

SET_TESTS_PROPERTIES(mrbiasAddMulMask PROPERTIES TIMEOUT 3600)
SET_TESTS_PROPERTIES(mrbiasMulIncremental PROPERTIES TIMEOUT 3600)
SET_TESTS_PROPERTIES(mrbiasMulLogIntensity PROPERTIES TIMEOUT 3600)

SET_TESTS_PROPERTIES(volume_reconstructionCubic PROPERTIES TIMEOUT 2400)
SET_TESTS_PROPERTIES(volume_reconstructionFourthOrder PROPERTIES TIMEOUT 1800)

SET_TESTS_PROPERTIES(warpDelayRefine PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(warpInverseConsistentCC PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(warpLabels PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(warpMultiLevel PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(warpMultiLevelExact PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(warpSingleLevel PROPERTIES TIMEOUT 1800)
SET_TESTS_PROPERTIES(warpSingleLevelExact PROPERTIES TIMEOUT 1800)

SET_TESTS_PROPERTIES(xform2scalarWarpInverseError PROPERTIES TIMEOUT 1800)
