#- Try to find the KinKal package
# Once done this will define
#  KINKAL_FOUND - System has KinKal
#  KINKAL_INCLUDE_DIR - The directory needed to compile against KinKal classes
#  KINKAL_LIBRARIES - The libraries needed to use KinKal


MESSAGE(STATUS "Adding KinKal include path $ENV{KINKAL_SOURCE_DIR}")
set(KINKAL_INCLUDE_DIR $ENV{KINKAL_SOURCE_DIR}/..)
mark_as_advanced(KINKAL_INCLUDE_DIR)
MESSAGE(STATUS "Looking for KinKal libraries in directory $ENV{KINKAL_LIBRARY_DIR}")
foreach( KKLIB MatEnv Trajectory Detector Fit General )
  set(KKFULLLIB KinKal_${KKLIB})
  MESSAGE("Looking for Kinkal library ${KKLIB}")
  find_library(${KKLIB} ${KKFULLLIB} HINTS $ENV{KINKAL_LIBRARY_DIR})
  if(${KKLIB})
    mark_as_advanced(${KKLIB})
    MESSAGE("Found KinKal library: ${KKLIB}")
    list(APPEND KINKAL_LIBRARIES ${KKFULLLIB})
  endif()
endforeach()

IF(KINKAL_LIBRARIES)
   SET(KINKAL_FOUND TRUE)
   MESSAGE(STATUS "Found KinKal libraries: ${KINKAL_LIBRARIES}")
ENDIF()

mark_as_advanced(KINKAL_LIBRARIES)

