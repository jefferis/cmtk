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

FIND_PROGRAM(SH_PATH sh)
IF(SH_PATH)
  SET(testDriver ${SH_PATH} ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.sh.in ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.sh @ONLY)
ELSE(SH_PATH)
  SET(testDriver python ${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.py)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/appsTestDriver.py.in ${CMAKE_CURRENT_BINARY_DIR}/appsTestDriver.py @ONLY)
ENDIF(SH_PATH)

# ==========================================
# Tests for "congeal" tool
SET(testList ${testList} 
  CongealFromInit
  CongealFromInitSampling
  CongealBackground 
  CongealUseTemplate
  CongealZeroSumSmooth)

# ==========================================
# Tests for "congeal_warp" tool
SET(testList ${testList} 
  CongealWarpFromInit
  CongealWarpFitFromInit
  CongealWarpFromInitZeroSum
  CongealWarpUseTemplate
  )

# ==========================================
# Tests for "describe" tool
SET(testList ${testList} DescribeMR1 DescribeMR2 DescribeMR3 DescribeMR4 DescribeMRBiorad DescribeNiftiDetached348)

IF(CMTK_BUILD_NRRD)
  SET(testList ${testList} DescribeMRNrrd1 DescribeMRNrrd2)
ENDIF(CMTK_BUILD_NRRD)

# ==========================================
# Tests for "groupwise_init" tool
SET(testList ${testList} GroupwiseInitCentersOfMass GroupwiseInitCentersOfMassTemplate)

# ==========================================
# Tests for "groupwise_rmi" tool
SET(testList ${testList}
  GroupwiseRMIFromInit
  GroupwiseRMIFromInitDeltaF
  GroupwiseRMIFromInitSampling 
  GroupwiseRMIBackground 
  GroupwiseRMIZeroSumSmooth)

IF(NOT CMTK_USE_MPI)
  # ==========================================
  # Tests for xml and wiki markups
  SET(testList_xml film levelset mrbias registration)
  SET(testList_wiki film levelset mrbias registration)

  # ==========================================
  # Tests for "concat_affine" tool
  SET(testList ${testList} 
    ConcatAffineABA
    ConcatAffineABAInvert
    ConcatAffineAB1A	     
    ConcatAffineAA1
    ConcatAffineA1A)

  # ==========================================
  # Tests for "average_edt" tool
  SET(testList ${testList} ShapeBasedAverage ShapeBasedAverageInterpolation)

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
    ConvertAutoCrop
    ConvertBoundaryMap
    ConvertBoundaryMapMultiValue
    ConvertDownsample
    ConvertErodeDilateErode
    ConvertDilateErodeDilate
    ConvertMedianFilter1
    ConvertMedianFilter2
    ConvertNiftiToAnalyze
    ConvertAnalyzeToNifti
    ConvertNiftiDetachedToNifti
    ConvertAnalyzeToNiftiDetached )

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      ConvertDownsampleNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "dcm2image" tool
  IF(CMTK_USE_DCMTK)
    SET(testList ${testList} Dcm2Image Dcm2ImageZ)

    IF(CMTK_BUILD_NRRD)
      SET(testList ${testList} Dcm2ImageNrrd)
    ENDIF(CMTK_BUILD_NRRD)
  ENDIF(CMTK_USE_DCMTK)

  # ==========================================
  # Tests for "filter" tool
  SET(testList ${testList} 
    FilterGaussian
    FilterGaussianSmallKernel
    FilterGaussianNoFilter)

  # ==========================================
  # Tests for "film" tool
  SET(testList ${testList} FilmCubic
    FilmFourthOrder
    FilmMSDLinearNoTrunc
    FilmMISincRefSPGR)

  # ==========================================
  # Tests for "gregxform" tool
  SET(testList ${testList} GregxformFordwardBackward GregxformAffine GregxformAffineFromWarp GregxformAffineFromWarpFwdBwd)

  # ==========================================
  # Tests for "histogram" tool
  SET(testList ${testList}
    Histogram 
    HistogramNorm
    HistogramBinsMinMax
    HistogramBinsMinMaxTrunc
    HistogramMask)

  # ==========================================
  # Tests for "imagemath" tool
  SET(testList ${testList} 
    ImagemathAverage ImagemathXor ImagemathContractLabels
    ImagemathStackEntropyLabels ImagemathVote ImagemathSTAPLE
    ImagemathMultiClassSTAPLE ImagemathCombinePCA ImagemathT2
    ImagemathLogOddsAdd ImagemathLogOddsAdd2
    ImagemathMatchHistograms ImagemathMatchHistogramsPadding)

  # ==========================================
  # Tests for "levelset" tool
  SET(testList ${testList} LevelsetDefault)

  # ==========================================
  # Tests for "make_initial_affine" tool
  SET(testList ${testList} 
    MakeInitialAffineCenterOfMass
    MakeInitialAffinePrincipalAxes1
    MakeInitialAffinePrincipalAxes2
    MakeInitialAffinePrincipalAxes3
    MakeInitialAffinePrincipalAxes4
    MakeInitialAffinePrincipalAxes5
    MakeInitialAffinePrincipalAxes6)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList}
      MakeInitialAffineDirectionVectorsNrrdAxSa
      MakeInitialAffineDirectionVectorsNrrdAxCo
      MakeInitialAffineDirectionVectorsNrrdSaCo
      MakeInitialAffineDirectionVectorsNrrdAxSaNative
      MakeInitialAffineDirectionVectorsNrrdAxCoNative
      MakeInitialAffineDirectionVectorsNrrdSaCoNative)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Multi-channel registration requires hash map to build
  IF(HAVE_STL_HASH_MAP)
    # ==========================================
    # Tests for "mcaffine" tool
    SET(testList ${testList} 
      McAffine1 McAffine2 McAffine3 McAffine4)
    
    # ==========================================
    # Tests for "mcwarp" tool
    SET(testList ${testList} 
      McWarp1 McWarp2 McWarp3)

  ENDIF(HAVE_STL_HASH_MAP)

  # ==========================================
  # Tests for "mk_phantom_3d" tool
  SET(testList ${testList} MkPhantomBox MkPhantomSphere MkPhantomBoxSphere MkPhantomImport)

  # ==========================================
  # Tests for "model" tool
  SET(testList ${testList} GLMDefault GLMNormalize GLMExp GLMNoConstant GLMIgnore GLMSelect GLMCrop)

  # ==========================================
  # Tests for "reorient" tool
  SET(testList ${testList}
    ReorientHdrSaToAx
    ReorientHdrCoToAx
    ReorientHdrAxToSa
    ReorientHdrCoToSa
    ReorientHdrAxToCo
    ReorientHdrSaToCo)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} ReorientNrrdToNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "split" tool
  SET(testList ${testList} SplitAxial SplitSagittal2 SplitCoronal3)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} SplitAxialNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "overlap" tool
  SET(testList ${testList}
    Overlap
    OverlapNumLabels
    OverlapByLabel
    OverlapFirst
    OverlapFirstByLabel
    OverlapFirstByLabelNumLabels)

  # ==========================================
  # Tests for "probe_xform" tool
  SET(testList ${testList}
    ProbeXformSimple
    ProbeXformFwdBwd
    ProbeXformBwdFwd)

  # ==========================================
  # Tests for "registration" tool
  SET(testList ${testList} 
    AffineRegistrationMrMrMSD
    RigidRegistrationMrPet 
    RigidRegistrationMrCt
    RigidRegistrationCt
    RigidRegistrationPetMr
    RigidRegistrationCtMr
    RigidRegistrationMrPetNoSwap
    RigidRegistrationMrCtNoSwap
    RigidRegistrationPetMrDOF9
    RigidRegistrationCtMrDOF7
    RigidRegistrationLabelsDOF69
    RigidRegistrationCrop
    RegistrationAutoLevelsRat4
    RegistrationAutoLevelsRat2
    RegistrationAutoLevelsRatToRat
    RegistrationAutoLevelsRatToRatDeltaFThreshold
    RegistrationAutoLevelsCt3
    RegistrationFromList
    RegistrationWithInitial
    RegistrationWithInitialInverse)

  # ==========================================
  # Tests for "mrbias" tool
  SET(testList ${testList} 
    MrbiasMulIncremental MrbiasMulLogIntensity MrbiasAddMulMask)

  # ==========================================
  # Tests for "reformatx" tool
  SET(testList ${testList} 
    ReformatxNoXform
    ReformatxLinear 
    ##		   ReformatxLinearFwdBwd
    ReformatxNearestNeighbor
    ReformatxPartialVolume
    ReformatxCubic
    ReformatxCubicInverse
    ReformatxSincCosine
    ReformatxSincCosine5
    ReformatxJacobian
    ReformatxInverseJacobian
    ReformatxSincHamming)
  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} ReformatxDfieldNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "sequence" tool
  SET(testList ${testList} SequenceDefault SequenceFormat SequenceThresh SequenceAbs SequenceAbsThresh)

  # ==========================================
  # Tests for "similarity" tool
  SET(testList ${testList} SimilarityGrey SimilarityLabels SimilarityGreyMask SimilarityLabelsMask)

  # ==========================================
  # Tests for "symmetry_plane" tool
  SET(testList ${testList} SymmetryPlane SymmetryPlaneThresh)

  # ==========================================
  # Tests for "statistics" tool
  SET(testList ${testList} 
    StatisticsGrey StatisticsPercentiles
    StatisticsGreyColumn StatisticsGreyExpNotation
    StatisticsGreyMask StatisticsGreyMultiMask	
    StatisticsLabels)

  # ==========================================
  # Tests for "ttest" tool
  SET(testList ${testList} TTestDefault TTestLog TTestAbsLog TTestInvert
    TTestPaired TTestCrossCorrelation TTestZScores
    TTestMask TTestOneSided)

  # ==========================================
  # Tests for "unsplit" tool
  SET(testList ${testList} UnsplitHdrAx UnsplitHdrSa UnsplitHdrCo)

  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} UnsplitHdrNrrdAx UnsplitHdrNrrdSa UnsplitHdrNrrdCo UnsplitNrrdNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "volume_injection" tool
  SET(testList ${testList} VolumeInjection VolumeInjectionIsotropic VolumeInjectionNoXform VolumeInjectionNoXformIsotropic)

  # ==========================================
  # Tests for "volume_reconstruction" tool
  SET(testList ${testList} VolumeReconstructionFourthOrder VolumeReconstructionCubic VolumeReconstructionNoXform)

  # ==========================================
  # Tests for "warp" tool
  SET(testList ${testList} 
    WarpSingleLevel
    WarpSingleLevelExact
    WarpInverseConsistentCC
    WarpMultiLevel
    WarpMultiLevelMatchHistograms
    WarpMultiLevelDeltaFThreshold
    WarpMultiLevelExact
    WarpDelayRefine
    WarpEnergy
    WarpJacobian
    WarpLabels)

  # ==========================================
  # Tests for "xform2dfield" tool
  IF(CMTK_BUILD_NRRD)
    SET(testList ${testList} 
      Xform2dfieldWarpNrrd
      Xform2dfieldAffineNrrd
      Xform2dfieldDownsampleXYZNrrd
      Xform2dfieldDownsampleXNrrd
      Xform2dfieldConcatNrrd
      Xform2dfieldInverseNrrd)
  ENDIF(CMTK_BUILD_NRRD)

  # ==========================================
  # Tests for "xform2scalar" tool
  SET(testList ${testList} 
    xform2scalarAffine 
    xform2scalarAffineDoubleY 
    xform2scalarWarp 
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

SET_TESTS_PROPERTIES(CongealWarpFromInit PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(CongealWarpFitFromInit PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(CongealWarpFromInitZeroSum PROPERTIES TIMEOUT 4800)
SET_TESTS_PROPERTIES(CongealWarpUseTemplate PROPERTIES TIMEOUT 6400)

IF(NOT CMTK_USE_MPI)
  # ==========================================
  # Set up test properties, now that the tests 
  # have been defined.

  SET_TESTS_PROPERTIES(FilmCubic PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(FilmMSDLinearNoTrunc PROPERTIES TIMEOUT 3600)
  SET_TESTS_PROPERTIES(FilmMISincRefSPGR PROPERTIES TIMEOUT 3600)

  SET_TESTS_PROPERTIES(McAffine2 PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(McWarp1 PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(McWarp2 PROPERTIES TIMEOUT 1800)

  SET_TESTS_PROPERTIES(MrbiasAddMulMask PROPERTIES TIMEOUT 3600)
  SET_TESTS_PROPERTIES(MrbiasMulIncremental PROPERTIES TIMEOUT 3600)
  SET_TESTS_PROPERTIES(MrbiasMulLogIntensity PROPERTIES TIMEOUT 3600)

  SET_TESTS_PROPERTIES(VolumeReconstructionCubic PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(VolumeReconstructionFourthOrder PROPERTIES TIMEOUT 1800)

  SET_TESTS_PROPERTIES(WarpDelayRefine PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(WarpInverseConsistentCC PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(WarpLabels PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(WarpMultiLevel PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(WarpMultiLevelExact PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(WarpSingleLevel PROPERTIES TIMEOUT 1800)
  SET_TESTS_PROPERTIES(WarpSingleLevelExact PROPERTIES TIMEOUT 1800)
ENDIF(NOT CMTK_USE_MPI)
