PROGRAMS = plot_decay cout_input

CXX=g++ 
CXXFLAGS= -Wno-deprecated -g -Wall -ggdb3 -fPIC -O2

CXXFLAGS += -I./include
ARFLAGS=


## for cern root
CXXFLAGS        += $(shell root-config --cflags)
LDFLAGS         += $(shell root-config --libs)
LDFLAGS         += $(shell root-config --glibs)

LIBNAME=lib/libTALYStool.a

AOBJS = ReadTALYS.o NucleusTable.o Nucleus.o

OBJDIR=obj
OBJS = $(addprefix $(OBJDIR)/,$(AOBJS))

EXEDIR=main
EXEOBJDIR=obj
EXEOBJS=$(addsuffix .o,$(addprefix $(EXEOBJDIR)/,$(PROGRAMS)))
BINDIR=bin
BIN=$(addprefix $(BINDIR)/,$(PROGRAMS))

all: $(BIN) lib

$(BIN): $(OBJS) $(EXEOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %$(notdir $@).o,$(EXEOBJS)) ${OBJS} $(LDFLAGS)

$(OBJDIR)/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXEOBJDIR)/%.o: main/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

lib:  $(OBJS)
	rm -rf $(LIBNAME)
	ar $(ARFLAGS) -rv $(LIBNAME) $^

clean: 
	rm -rf $(OBJDIR)/*.o $(BINDIR)/* $(LIBNAME)
