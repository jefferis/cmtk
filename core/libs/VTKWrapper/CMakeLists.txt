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

INSTALL(TARGETS cmtkVTKWrapper
  RUNTIME DESTINATION bin COMPONENT runtime
  LIBRARY DESTINATION lib COMPONENT runtime
  ARCHIVE DESTINATION lib COMPONENT libraries)

FILE(GLOB files_h "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
FILE(GLOB files_txx "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${files_h} ${files_txx} DESTINATION include/cmtk COMPONENT development)
