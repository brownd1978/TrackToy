# generate root dictionary
ROOT_GENERATE_DICTIONARY(TestDict
    # headers to include in ROOT dictionary
    TrkInfo.hh
    LINKDEF LinkDef.h
    MODULE Test
)
# create shared library with ROOT dict
add_library(Test SHARED TestDict)
target_include_directories(Test PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# link ROOT libraries
target_link_libraries(Test ${ROOT_LIBRARIES})
