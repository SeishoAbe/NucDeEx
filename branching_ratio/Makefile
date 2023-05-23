PROGRAMS = plot_decay

CXX=g++ 
CXXFLAGS= -Wno-deprecated -g -Wall -ggdb3 -fPIC -O2

CXXFLAGS += -I./include
ARFLAGS=


## for cern root
CXXFLAGS        += $(shell root-config --cflags)
LDFLAGS         += $(shell root-config --libs)
LDFLAGS         += $(shell root-config --glibs)

LIBNAME=lib/libNucleusTable.a

AOBJS = NucleusTable.o

OBJDIR=obj
OBJS = $(addprefix $(OBJDIR)/,$(AOBJS))

EXEDIR=main
EXEOBJDIR=main
EXEOBJ=$(EXEOBJDIR)/$(PROGRAMS).o
BINDIR=bin
BIN=$(addprefix $(BINDIR)/,$(PROGRAMS))

all: $(BIN) lib

$(BIN): $(OBJS) $(EXEOBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(EXEOBJ) ${OBJS} $(LDFLAGS)

$(OBJDIR)/%.o: src/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXEOBJDIR)/%.o: main/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

lib:  $(OBJS)
	rm -rf $(LIBNAME)
	ar $(ARFLAGS) -rv $(LIBNAME) $^

clean: 
	rm -rf $(OBJDIR)/*.o $(BINDIR)/* $(LIBNAME)
