# Sources of non-templated classes.
SET(cmtkVTKWrapper_SRCS
  cmtkLandmarkListVTK.cxx
  cmtkMatchedLandmarkListVTK.cxx
  cmtkVolumeToVtkStructuredPoints.cxx
  cmtkVtkADM.cxx
  cmtkVtkIsoSurface.cxx
  cmtkVtkLookupTable.cxx
  cmtkWarpVtkPolyData.cxx
)

ADD_LIBRARY(cmtkVTKWrapper ${cmtkVTKWrapper_SRCS})

TARGET_LINK_LIBRARIES(cmtkVTKWrapper cmtkIO cmtkBase cmtkSystem cmtkNumerics )

IF(CMTK_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(cmtkVTKWrapper PROPERTIES ${CMTK_LIBRARY_PROPERTIES})
ENDIF(CMTK_LIBRARY_PROPERTIES)

IF(${CMTK_BUILD_SHARED_LIBS})
  INSTALL(TARGETS cmtkVTKWrapper LIBRARY DESTINATION bin)
ENDIF(${CMTK_BUILD_SHARED_LIBS})
