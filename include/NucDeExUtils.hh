#ifndef __NUCDEEXUTILS__HH__
#define __NUCDEEXUTILS__HH__

#include <TDatabasePDG.h>
#include <TGeoManager.h>
#include <TGeoElement.h>

#include "NucDeExNucleusTable.hh"

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

namespace NucDeEx{
  namespace Utils{
    extern int fVerbose;
    extern NucDeExNucleusTable* NucleusTable;
    extern std::string NUCDEEX_ROOT;
    extern TDatabasePDG* fTDatabasePDG;
    extern TGeoManager* fTGeoManager;
    extern TGeoElementTable* fTGeoElementTable;
    // needs "extern". These are defined in *.cc

    void Init();

    void SetPATH();
#ifdef INCL_DEEXCITATION_NUCDEEX
    void SetPATH(G4INCL::Config* config);
#endif
  }
}
#endif
