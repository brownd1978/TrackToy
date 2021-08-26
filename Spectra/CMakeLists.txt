
add_library(Spectra SHARED 
    CeSpectrum.cc
)

# set top-level directory as include root
target_include_directories(Spectra PRIVATE ${PROJECT_SOURCE_DIR}/..)

# set shared library version equal to project version
set_target_properties(Spectra PROPERTIES VERSION ${PROJECT_VERSION})
