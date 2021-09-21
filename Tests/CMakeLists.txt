# List of test test sources
# ls -1 *_test.cc

set( TEST_SOURCE_FILES
  CeMinusSpectrum_test.cc
  DIOSpectrum_test.cc
  MuonRange_test.cc
  Particle_test.cc
  BField_test.cc
  MuStops_test.cc
)

# generate root dictionary
include_directories(${PROJECT_SOURCE_DIR}/..)
ROOT_GENERATE_DICTIONARY(TestsDict
    ${PROJECT_SOURCE_DIR}/../KinKal/General/ParticleState.hh
    LINKDEF LinkDef.h
)
# create shared library with ROOT dict
add_library(Tests SHARED TestsDict)

target_include_directories(Tests PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# link ROOT libraries
target_link_libraries(Tests ${ROOT_LIBRARIES} ${KINKAL_LIBRARIES})

# Generate test targets
foreach( testsourcefile ${TEST_SOURCE_FILES} )
    # get the name of the test test from the source path
    # (strip _test.cc, and the path preceding the basename)
    string( REPLACE "_test.cc" "" testnamenoext ${testsourcefile} )
    get_filename_component(testname ${testnamenoext} NAME)

    # prepend Test_ to the target name to avoid target name clashes,
    # create the test test executable:
    add_executable( Test_${testname} ${testsourcefile} )

    # add the project root as an include directory
    # link all test tests to KinKal, MatEnv, and ROOT libraries.
    target_include_directories(Test_${testname} PRIVATE ${PROJECT_SOURCE_DIR}/..)
    target_link_libraries( Test_${testname}  Spectra Detector General ${ROOT_LIBRARIES} ${KINKAL_LIBRARIES})

    # ensure the test test executable filename is just its test name
    set_target_properties( Test_${testname} PROPERTIES OUTPUT_NAME ${testname})

    # register the target as a test
    add_test (NAME ${testname} COMMAND Test_${testname} )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 250)
    set_tests_properties(${testname} PROPERTIES ENVIRONMENT "TRACKTOY_SOURCE_DIR=${CMAKE_SOURCE_DIR}")

    install( TARGETS Test_${testname}
             RUNTIME DESTINATION bin/ )

endforeach( testsourcefile ${TEST_SOURCE_FILES} )
