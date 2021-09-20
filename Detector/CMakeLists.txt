
add_library(Detector SHARED
#    EdmlXMLEngine.cc
#   EdmlParser.cc
)

# set top-level directory as include root
target_include_directories(Detector PRIVATE ${PROJECT_SOURCE_DIR}/..)

# set shared library version equal to project version
set_target_properties(Detector PROPERTIES VERSION ${PROJECT_VERSION})
