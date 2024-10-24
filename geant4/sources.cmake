#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_hetcpp_utils
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_hetcpp.G4hadronic_hetcpp_utils
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List root includes needed.
execute_process(
  COMMAND root-config --incdir
  OUTPUT_VARIABLE ROOT_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
include_directories(${ROOT_INCLUDE_DIR})

execute_process(
  COMMAND root-config --libs
  OUTPUT_VARIABLE ROOT_LIBRARIES
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lEG -lGeom) # -lEve will break G4 if it's linked
#message(${ROOT_LIBRARIES})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/handler/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_nucdeex
    HEADERS
        G4NucDeExInterface.hh
        G4NucDeExInterfaceMessenger.hh
        NucDeExConsts.hh
        NucDeExDeexcitation.hh
        NucDeExDeexcitationBase.hh
        NucDeExDeexcitationPhole.hh
        NucDeExDeexcitationTALYS.hh
        NucDeExEventInfo.hh
        NucDeExNucleus.hh
        NucDeExNucleusTable.hh
        NucDeExParticle.hh
        NucDeExRandom.hh
        NucDeExUtils.hh

    SOURCES
        G4NucDeExInterface.cc
        G4NucDeExInterfaceMessenger.cc
        NucDeExDeexcitation.cc
        NucDeExDeexcitationBase.cc
        NucDeExDeexcitationPhole.cc
        NucDeExDeexcitationTALYS.cc
        NucDeExEventInfo.cc
        NucDeExNucleus.cc
        NucDeExNucleusTable.cc
        NucDeExParticle.cc
        NucDeExRandom.cc
        NucDeExUtils.cc

    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4hadronic_mgt
        G4hadronic_util
        G4hadronic_xsect
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4track
        G4volumes
        G4intercoms
        G4had_preequ_exciton
        G4had_mod_man
        G4had_mod_util
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util

    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
        G4intercoms

    LINK_LIBRARIES
      ${ROOT_LIBRARIES}
)# List any source specific properties here
