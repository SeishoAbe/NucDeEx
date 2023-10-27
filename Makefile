PROGRAMS = plot_decay cout_input simulation
PROGRAMS += genie

CXX=g++ 
CXXFLAGS= -Wno-deprecated -g -Wall -ggdb3 -fPIC -O2

CXXFLAGS += -I./include
ARFLAGS=

## for cern root
CXXFLAGS        += $(shell root-config --cflags)
#LDFLAGS         += $(shell root-config --libs) # <- Print regular ROOT libraries
#LDFLAGS         += $(shell root-config --glibs) # <- Print regular + GUI ROOT libraries
LDFLAGS         += $(shell root-config --evelibs) # <- Print regular + GUI + Eve libraries. SHOULD BE THIS!

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
LIBNAME=${LIBDIR}/libTALYStool.a

AOBJS =  Nucleus.o NucleusTable.o ReadTALYS.o 
AOBJS += Particle.o Deexcitation.o

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
