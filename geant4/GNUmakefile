# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4hadronic_nucdeex

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/geometry/management/include \
	    -I$(G4BASE)/geometry/volumes/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/materials/include \
	    -I$(G4BASE)/particles/bosons/include \
	    -I$(G4BASE)/particles/hadrons/barions/include \
	    -I$(G4BASE)/particles/hadrons/ions/include \
	    -I$(G4BASE)/particles/hadrons/mesons/include \
	    -I$(G4BASE)/particles/leptons/include \
	    -I$(G4BASE)/particles/management/include \
	    -I$(G4BASE)/processes/management/include \
	    -I$(G4BASE)/processes/hadronic/models/management/include \
	    -I$(G4BASE)/processes/hadronic/util/include \
	    -I$(G4BASE)/processes/hadronic/models/util/include \
	    -I$(G4BASE)/processes/hadronic/models/de_excitation/handler/include \
	    -I$(G4BASE)/processes/hadronic/models/de_excitation/evaporation/include \
	    -I$(G4BASE)/track/include \
	    -I$(G4BASE)/intercoms/include
           
include $(G4INSTALL)/config/common.gmk
