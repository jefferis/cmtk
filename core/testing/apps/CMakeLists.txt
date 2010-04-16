##
##  Copyright 2004-2010 SRI International
##
##  Copyright 1997-2009 Torsten Rohlfing
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

FIND_PROGRAM(SH_PATH sh)
IF(SH_PATH)
  SET(testDriver ${SH_PATH} ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.sh.in ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh @ONLY)
ELSE(SH_PATH)
  SET(testDriver python ${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.py)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.py.in ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.py @ONLY)
ENDIF(SH_PATH)


# ==========================================
# Tests for "avg_adm" tool
SET(testList ${testList} 
  avg_admDefault2
  avg_admDefault3
  avg_admJacobianFloat
  avg_admLabels
  avg_admPaddOutZeroNN
  avg_admNoReferenceModelCubic
  avg_admNoScaleModelAutoScale
  avg_admWithAffineNoRefData)

# ==========================================
# Tests for "congeal" tool
SET(testList ${testList} 
  congealFromInit
  congealFromInitSampling
  congealBackground 
  congealUseTemplate
  congealZeroSumSmooth)

# ==========================================
# Tests for "congeal_warp" tool
SET(testList ${testList} 
  congeal_warpFromInit
  congeal_warpFitFromInit
  congeal_warpFromInitZeroSum
  congeal_warpUseTemplate
  )

# ==========================================
# Tests for "dbtool"
IF(CMTK_BUILD_SQLITE)
  SET(testList ${testList}
    dbtool_AddImages
    dbtool_AddImages2
    dbtool_ListSpace
    dbtool_GetXform1
    dbtool_GetXform2
    dbtool_GetXform3
    dbtool_GetXform4
    dbtool_GetXform5
    dbtool_GetXform6)
ENDIF(CMTK_BUILD_SQLITE)  
  
# ==========================================
# Tests for "describe" tool
SET(testList ${testList} 
  describeEmpty
  describeMR1 
  describeMR2
  describeMR3
  describeMR4
  describeMRBiorad
  describeNiftiDetached348)

IF(CMTK_BUILD_NRRD)
  SET(testList ${testList} describeMRNrrd1 describeMRNrrd2)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "groupwise_init" tool
SET(testList ${testList} groupwise_initCentersOfMass groupwise_initCentersOfMassTemplate)

# ==========================================
# Tests for "groupwise_rmi" tool
SET(testList ${testList}
  groupwise_rmiFromInit
  groupwise_rmiFromInitDeltaF
  groupwise_rmiFromInitSampling 
  groupwise_rmiBackground 
  groupwise_rmiZeroSumSmooth)

IF(NOT CMTK_USE_MPI)
  # ==========================================
  # Tests for xml and wiki markups
  SET(testList_xml film levelset mrbias registration)
  SET(testList_wiki film levelset mrbias registration)

  # ==========================================
  # Tests for "concat_affine" tool
  SET(testList ${testList} 
    concat_affineABA
    concat_affineABAInvert
    concat_affineAB1A	     
    concat_affineAA1
    concat_affineA1A)

  # ==========================================
  # Tests for "average_edt" tool
  SET(testList ${testList} 
    average_edt 
    average_edtWindowed
    average_edtWindowedInterpolation
    average_edtInterpolation)

  # ==========================================
  # Tests for "average_images" tool
  SET(testList ${testList} 
    average_imagesMean 
    average_imagesMeanNormPadd 
    average_imagesMeanAbsLog 
    average_imagesVariance
    average_imagesStDev
    average_imagesEntropy
    average_imagesZScore)

  # ==========================================
  # Tests for "convert" tool
  SET(testList ${testList} 
    convertAutoCrop
    convertBoundaryMap
    convertBoundaryMapMultiValue
    convertDownsample
    convertErodeDilateErode
    convertDilateErodeDilate
    convertMedianFilter1
    convertMedianFilter2
    convertNiftiToAnalyze
    convertAnalyzeToNifti
    convertNiftiDetachedToNifti
    convertAnalyzeToNiftiDetached )

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      convertDownsampleNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "convertx" tool
  SET(testList ${testList} 
    convertxBinarize
    convertxBoundaryMap
    convertxBoundaryMapMultiValue
    convertxCropThresholdWriteRegion
    convertxCropRegion1
    convertxCropRegion2
    convertxDownsample
    convertxDistanceUnsigned
    convertxDistanceSigned
    convertxDistanceSigned2
    convertxErodeDilateErode
    convertxDilateErodeDilate
    convertxGaussianFilterSigma
    convertxGaussianFilterFWHM
    convertxHistogramEqualization
    convertxHistogramEqualizationNBins
    convertxMask
    convertxMaskInverse
    convertxMedianFilter1
    convertxMedianFilter2
    convertxMedianFilterXYZ
    convertxNiftiToAnalyze
    convertxAnalyzeToNifti
    convertxNiftiDetachedToNifti
    convertxAnalyzeToNiftiDetached 
    convertxThresholdBelow
    convertxThresholdAbove
    convertxThresholdBelowToPadding
    convertxThresholdAboveToPadding
    convertxScaleToRange
    convertxRevert)
  
  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      convertxDownsampleNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "dcm2image" tool
  IF(CMTK_USE_DCMTK)
    SET(testList ${testList} dcm2image dcm2imageZ)

    IF(CMTK_BUILD_NRRD)
      SET(testList ${testList} dcm2imageNrrd)
    ENDIF(CMTK_BUILD_NRRD)
  ENDIF(CMTK_USE_DCMTK)

  # ==========================================
  # Tests for "filter" tool
  SET(testList ${testList} 
    filterGaussian
    filterGaussianSmallKernel
    filterGaussianNoFilter
    filterRohlfing)

  # ==========================================
  # Tests for "film" tool
  SET(testList ${testList} filmCubic
    filmFourthOrder
    filmMSDLinearNoTrunc
    filmMISincRefSPGR)

  # ==========================================
  # Tests for "gregxform" tool
  SET(testList ${testList} gregxformFordwardBackward gregxformAffine gregxformAffineFromWarp gregxformAffineFromWarpFwdBwd)

  # ==========================================
  # Tests for "histogram" tool
  SET(testList ${testList}
    histogram 
    histogramNorm
    histogramBinsMinMax
    histogramBinsMinMaxTrunc
    histogramMask)

  # ==========================================
  # Tests for "imagemath" tool
  SET(testList ${testList} 
    imagemathAverage 
    imagemathXor
    imagemathDupFillMaxValue
    imagemathAnd
    imagemathContractLabels
    imagemathStackEntropyLabels
    imagemathVote
    imagemathSTAPLE
    imagemathMultiClassSTAPLE 
    imagemathCombinePCA
    imagemathT2
    imagemathLogOddsAdd
    imagemathLogOddsAdd2
    imagemathMatchHistograms 
    imagemathMatchHistogramsPadding
    imagemathMatchHistogramsPadding2
    imagemathMatchHistogramsPaddingUnset)

  # ==========================================
  # Tests for "levelset" tool
  SET(testList ${testList} levelsetDefault)

  # ==========================================
  # Tests for "make_initial_affine" tool
  SET(testList ${testList} 
    make_initial_affineCenterOfMass
    make_initial_affinePrincipalAxes1
    make_initial_affinePrincipalAxes2
    make_initial_affinePrincipalAxes3
    make_initial_affinePrincipalAxes4
    make_initial_affinePrincipalAxes5
    make_initial_affinePrincipalAxes6)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList}
      make_initial_affineDirectionVectorsNrrdAxSa
      make_initial_affineDirectionVectorsNrrdAxCo
      make_initial_affineDirectionVectorsNrrdSaCo
      make_initial_affineDirectionVectorsNrrdAxSaNative
      make_initial_affineDirectionVectorsNrrdAxCoNative
      make_initial_affineDirectionVectorsNrrdSaCoNative)
  ENDIF(CMTK_BUILD_NRRD)

  IF(CMTK_BUILD_SQLITE)
    SET(testList ${testList}
      make_initial_affineDB)
  ENDIF(CMTK_BUILD_SQLITE)  
  
  # ==========================================
  # Tests for "mat2dof" tool
  SET(testList ${testList} 
    mat2dof1
    mat2dof2
    mat2dofFile)

  # ==========================================
  # Multi-channel registration requires hash map to build
  IF(HAVE_STL_HASH_MAP)
    # ==========================================
    # Tests for "mcaffine" tool
    SET(testList ${testList} 
      mcaffine1 mcaffine2 mcaffine3 mcaffine4)
    
    # ==========================================
    # Tests for "mcwarp" tool
    SET(testList ${testList} 
      mcwarp1 mcwarp2 mcwarp3)

  ENDIF(HAVE_STL_HASH_MAP)

  # ==========================================
  # Tests for "mip" tool
  SET(testList ${testList} 
    mipDefault
    mipBlackWhite
    mipAdjust
    mipSum
    mip16Bit)

  # ==========================================
  # Tests for "mk_analyze_hdr" tool
  SET(testList ${testList} 
    mk_analyze_hdrDefault
    mk_analyze_hdrImport
    mk_analyze_hdrLittleEndian
    mk_analyze_hdrBigEndian)

  # ==========================================
  # Tests for "mk_phantom_3d" tool
  SET(testList ${testList} 
    mk_phantom_3dBoxIndexed 
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
  # Tests for "model" tool
  SET(testList ${testList} glmDefault glmNormalize glmExp glmNoConstant glmIgnore glmSelect glmCrop)

  # ==========================================
  # Tests for "reorient" tool
  SET(testList ${testList}
    reorientHdrSaToAx
    reorientHdrCoToAx
    reorientHdrAxToSa
    reorientHdrCoToSa
    reorientHdrAxToCo
    reorientHdrSaToCo)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} reorientNrrdToNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "split" tool
  SET(testList ${testList} splitAxial splitSagittal2 splitCoronal3)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} splitAxialNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "overlap" tool
  SET(testList ${testList}
    overlap
    overlapNumLabels
    overlapByLabel
    overlapFirst
    overlapFirstByLabel
    overlapFirstByLabelNumLabels)

  # ==========================================
  # Tests for "probe_xform" tool
  SET(testList ${testList}
    probe_xformSimple
    probe_xformFwdBwd
    probe_xformBwdFwd)

  # ==========================================
  # Tests for "registration" tool
  SET(testList ${testList} 
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
  SET(testList ${testList} 
    registrationxAffineMrMrMSD
    registrationxFromList
    registrationxWithInitial
    registrationxWithInitialInverse
    registrationxAutoLevelsRat4
    registrationxAutoLevelsRat4COM
    registrationxAutoLevelsRat4PAX
    registrationxAutoLevelsRat2
    registrationxAutoLevelsRat2Cubic
    registrationxAutoLevelsRat2Sinc
    registrationxAutoLevelsRat2NN
    registrationxAutoLevelsRatToRat
    registrationxAutoLevelsRatToRatDeltaFThreshold
    registrationxAutoLevelsRatToRatNCC
    registrationxAutoLevelsRatToRatCR
    registrationxAutoLevelsRatToRatMI
    registrationxAutoLevelsLabelsNN
    registrationxAutoLevelsLabelsPV
    registrationxAutoLevelsCt3)

  IF(CMTK_BUILD_SQLITE)
    SET(testList ${testList} 
      registrationxFromListDB
      registrationxWithInitialDB
      registrationxWithInitialInverseDB
      )
  ENDIF(CMTK_BUILD_SQLITE)  
  
  # ==========================================
  # Tests for "mrbias" tool
  SET(testList ${testList} 
    mrbiasMulIncremental mrbiasMulLogIntensity mrbiasAddMulMask)

  # ==========================================
  # Tests for "reformatx" tool
  SET(testList ${testList} 
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
    reformatxTargetGridOffset)
  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      reformatxDfieldNrrd
      reformatxDfieldNrrdJacobian)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "sequence" tool
  SET(testList ${testList} sequenceDefault sequenceFormat sequenceThresh sequenceAbs sequenceAbsThresh)

  # ==========================================
  # Tests for "similarity" tool
  SET(testList ${testList} similarityGrey similarityLabels similarityGreyMask similarityLabelsMask)

  # ==========================================
  # Tests for "symmetry_plane" tool
  SET(testList ${testList} symmetry_plane symmetry_planeThresh)

  # ==========================================
  # Tests for "statistics" tool
  SET(testList ${testList} 
    statisticsGrey 
    statisticsPercentiles
    statisticsGreyColumn 
    statisticsGreyExpNotation
    statisticsGreyMask 
    statisticsGreyMultiMask
    statisticsMaskMismatch
    statisticsLabels)

  # ==========================================
  # Tests for "ttest" tool
  SET(testList ${testList} ttestDefault ttestLog ttestAbsLog ttestInvert
    ttestPaired ttestCrossCorrelation ttestZScores
    ttestMask ttestOneSided)

  # ==========================================
  # Tests for "unsplit" tool
  SET(testList ${testList} unsplitHdrAx unsplitHdrSa unsplitHdrCo)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} unsplitHdrNrrdAx unsplitHdrNrrdSa unsplitHdrNrrdCo unsplitNrrdNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "volume_injection" tool
  SET(testList ${testList} volume_injection volume_injectionIsotropic volume_injectionNoXform volume_injectionNoXformIsotropic)

  # ==========================================
  # Tests for "volume_reconstruction" tool
  SET(testList ${testList} volume_reconstructionFourthOrder volume_reconstructionCubic volume_reconstructionNoXform)

  # ==========================================
  # Tests for "warp" tool
  SET(testList ${testList} 
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
  SET(testList ${testList} 
    warpxLabels)

  # ==========================================
  # Tests for "xform2dfield" tool
  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      xform2dfieldWarpNrrd
      xform2dfieldAffineNrrd
      xform2dfieldDownsampleXYZNrrd
      xform2dfieldDownsampleXNrrd
      xform2dfieldConcatNrrd
      xform2dfieldInverseNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "xform2scalar" tool
  SET(testList ${testList} 
    xform2scalarAffine 
    xform2scalarAffineDoubleY 
    xform2scalarWarp 
    xform2scalarWarpInverseError
    xform2scalarWarpOnly)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      xform2scalarDfield)
  ENDIF(CMTK_BUILD_NRRD)
  
ENDIF(NOT CMTK_USE_MPI)

# ==========================================
# Set up all tests
FOREACH(testName ${testList})
  IF(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(${testName} ${testDriver} ${testName} ${MEMORYCHECK_COMMAND})
  ELSE(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(${testName} ${testDriver} ${testName})
  ENDIF(CMTK_TESTING_MEMORYCHECK)
ENDFOREACH(testName ${testList})

FOREACH(testName ${testList_xml})
  IF(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(xml_${testName} ${testDriver} xml_${testName} ${MEMORYCHECK_COMMAND})
  ELSE(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(xml_${testName} ${testDriver} xml_${testName})
  ENDIF(CMTK_TESTING_MEMORYCHECK)
ENDFOREACH(testName ${testList_xml})

FOREACH(testName ${testList_wiki})
  IF(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(wiki_${testName} ${testDriver} wiki_${testName} ${MEMORYCHECK_COMMAND})
  ELSE(CMTK_TESTING_MEMORYCHECK)
    ADD_TEST(wiki_${testName} ${testDriver} wiki_${testName})
  ENDIF(CMTK_TESTING_MEMORYCHECK)
ENDFOREACH(testName ${testList_wiki})

SET_TESTS_PROPERTIES(congeal_warpFromInit PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(congeal_warpFitFromInit PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(congeal_warpFromInitZeroSum PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(congeal_warpUseTemplate PROPERTIES TIMEOUT 6400)

IF(NOT CMTK_USE_MPI)
  # ==========================================
  # Set up test properties, now that the tests 
  # have been defined.

  SET_TESTS_PROPERTIES(filmCubic PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(filmMSDLinearNoTrunc PROPERTIES TIMEOUT 3600)
  SET_TESTS_PROPERTIES(filmMISincRefSPGR PROPERTIES TIMEOUT 3600)

  IF(HAVE_STL_HASH_MAP)
    SET_TESTS_PROPERTIES(mcaffine2 PROPERTIES TIMEOUT 1800)
    SET_TESTS_PROPERTIES(mcwarp1 PROPERTIES TIMEOUT 1800)
    SET_TESTS_PROPERTIES(mcwarp2 PROPERTIES TIMEOUT 1800)
  ENDIF(HAVE_STL_HASH_MAP)

  SET_TESTS_PROPERTIES(mrbiasAddMulMask PROPERTIES TIMEOUT 3600)
  SET_TESTS_PROPERTIES(mrbiasMulIncremental PROPERTIES TIMEOUT 3600)
  SET_TESTS_PROPERTIES(mrbiasMulLogIntensity PROPERTIES TIMEOUT 3600)

  SET_TESTS_PROPERTIES(volume_reconstructionCubic PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(volume_reconstructionFourthOrder PROPERTIES TIMEOUT 1800)

  SET_TESTS_PROPERTIES(warpDelayRefine PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(warpInverseConsistentCC PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(warpLabels PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(warpMultiLevel PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(warpMultiLevelExact PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(warpSingleLevel PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(warpSingleLevelExact PROPERTIES TIMEOUT 1800)

  SET_TESTS_PROPERTIES(xform2scalarWarpInverseError PROPERTIES TIMEOUT 1800)
ENDIF(NOT CMTK_USE_MPI)
