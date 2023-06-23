PROGRAMS = plot_decay cout_input simulation

CXX=g++ 
CXXFLAGS= -Wno-deprecated -g -Wall -ggdb3 -fPIC -O2

CXXFLAGS += -I./include
ARFLAGS=

## for cern root
CXXFLAGS        += $(shell root-config --cflags)
#LDFLAGS         += $(shell root-config --libs) # <- Print regular ROOT libraries
#LDFLAGS         += $(shell root-config --glibs) # <- Print regular + GUI ROOT libraries
LDFLAGS         += $(shell root-config --evelibs) # <- Print regular + GUI + Eve libraries. SHOULD BE THIS!

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
	rm -rf $(LIBNAME)
	ar $(ARFLAGS) -rv $(LIBNAME) $^

dylib:  $(OBJS)
	rm -rf $(LIBNAME:.a=.so)
	$(CXX) -shared -o $(LIBNAME:.a=.so) $^

clean: 
	rm -rf $(OBJDIR)/*.o $(BINDIR)/* $(LIBNAME)
