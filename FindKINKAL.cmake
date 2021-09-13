#- Try to find the KinKal package
# Once done this will define
#  KINKAL_FOUND - System has KinKal
#  KINKAL_LIBRARIES - The libraries needed to use KinKal

MESSAGE(STATUS "Looking for KinKal libraries in directory $ENV{KINKAL_LIBRARY_DIR}")
foreach( _cpt General MatEnv Trajectory)
  MESSAGE("Looking for Kinkal library ${_cpt}")
  find_library(KINKAL_${_cpt}_LIBRARY ${_cpt} HINTS $ENV{KINKAL_LIBRARY_DIR})
  if(KINKAL_${_cpt}_LIBRARY)
    mark_as_advanced(KINKAL_${_cpt}_LIBRARY)
    MESSAGE("Found KinKal library: ${KINKAL_${_cpt}_LIBRARY}")
    list(APPEND KINKAL_LIBRARIES ${KINKAL_${_cpt}_LIBRARY})
  endif()
endforeach()

IF(KINKAL_LIBRARIES)
   SET(KINKAL_FOUND TRUE)
   MESSAGE(STATUS "Found KinKal libraries: ${KINKAL_LIBRARIES}")
ENDIF()

mark_as_advanced(KINKAL_LIBRARIES)

