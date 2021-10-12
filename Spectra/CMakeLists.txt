add_library(Spectra SHARED
    CeMinusSpectrum.cc
    DIOSpectrum.cc
)

# set top-level directory as include root
target_include_directories(Spectra PRIVATE ${PROJECT_SOURCE_DIR}/..)
# add KinKal; this may be redundant, but...
target_include_directories(Spectra PRIVATE ${KINKAL_INCLUDE_DIR})

# set shared library version equal to project version
set_target_properties(Spectra PROPERTIES VERSION ${PROJECT_VERSION})
