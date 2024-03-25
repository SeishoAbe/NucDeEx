#ifndef __NUCDEEXUTILS__HH__
#define __NUCDEEXUTILS__HH__

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif
#include "NucDeExNucleusTable.hh"

namespace NucDeEx{
  namespace Utils{
    extern int fVerbose;
    extern NucDeExNucleusTable* NucleusTable;
    extern std::string NUCDEEX_ROOT;
    // needs "extern". These are defined in *.cc

    void SetPATH();
#ifdef INCL_DEEXCITATION_NUCDEEX
    void SetPATH(G4INCL::Config* config);
#endif
  }
}
#endif
