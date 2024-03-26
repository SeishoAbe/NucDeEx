PROGRAMS = plot_decay cout_input simulation
PROGRAMS += genie

ifndef ROOTSYS
@echo "$(error "Please set ROOT enviromental variables")"
endif

CXX=g++ 
CXXFLAGS= -Wno-deprecated -g -Wall -ggdb3 -fPIC -O2

CXXFLAGS += -I./include
ARFLAGS=

## for cern root
CXXFLAGS        += $(shell root-config --cflags)
LDFLAGS         += $(shell root-config --libs) -lEG -lGeom #-lEve 
# DO NOT link lEve. -lEG and -lGeom are necessary.
ROOTVERSION = $(shell root-config --version)
ROOT5=$(findstring 5.,$(ROOTVERSION))
ifneq ($(ROOT5),)
CXXFLAGS += -DROOT5
endif

## neut libs
ifdef NEUT_ROOT
PROGRAMS += neut
CXXFLAGS += -I$(NEUT_ROOT)/include
LDFLAGS += -L$(NEUT_ROOT)/lib -lNEUT -lNEUTClass -lNEUTClassUtils
endif

## nuwro libs
ifdef NUWRO
PROGRAMS += nuwro
CXXFLAGS += -I$(NUWRO)/src
LDFLAGS += $(NUWRO)/bin/event1.so
endif

LIBDIR=lib
LIBNAME=${LIBDIR}/libNucDeEx.a

AOBJS = ReadTALYS.o
AOBJS += NucDeExNucleus.o NucDeExNucleusTable.o
AOBJS += NucDeExParticle.o NucDeExUtils.o NucDeExEventInfo.o NucDeExRandom.o
AOBJS += NucDeExDeexcitationBase.o NucDeExDeexcitationTALYS.o NucDeExDeexcitationPhole.o NucDeExDeexcitation.o

OBJDIR=obj
OBJS = $(addprefix $(OBJDIR)/,$(AOBJS))

EXEDIR=main
EXEOBJDIR=obj
EXEOBJS=$(addsuffix .o,$(addprefix $(EXEOBJDIR)/,$(PROGRAMS)))
BINDIR=bin
BIN=$(addprefix $(BINDIR)/,$(PROGRAMS))

all: $(BIN) lib dylib

$(BIN): $(OBJS) $(EXEOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %$(notdir $@).o,$(EXEOBJS)) ${OBJS} $(LDFLAGS)

$(OBJDIR)/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXEOBJDIR)/%.o: main/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

lib:  $(OBJS)
	ar $(ARFLAGS) -rv $(LIBNAME) $^

dylib:  $(OBJS)
	$(CXX) -shared -o $(LIBNAME:.a=.so) $^

clean: 
	rm -rf $(OBJDIR)/*.o $(BINDIR)/* $(LIBDIR)/*
