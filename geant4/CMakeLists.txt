#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4hadronic_nucdeex
# Package:
# Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_nucdeex
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

