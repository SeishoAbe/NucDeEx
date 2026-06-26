#include <iostream>
#include <string>
#include <cstdlib>

#include "NucDeExUtils.hh"

namespace NucDeEx{
  namespace Utils{
    // --- Define & initialize --- //
    int fVerbose=0;
    std::string NUCDEEX_ROOT="";
    NucDeExNucleusTable* NucleusTable=nullptr;
    TDatabasePDG* fTDatabasePDG=nullptr;
    TGeoManager* fTGeoManager=nullptr;
    TGeoElementTable* fTGeoElementTable=nullptr;

    void Init()
    {
      if(NucleusTable==nullptr){
         NucleusTable = new NucDeExNucleusTable();
      }
      if(fTDatabasePDG==nullptr){
        fTDatabasePDG = new TDatabasePDG();
      }
      std::cout << "gGeoManager=" << gGeoManager << std::endl;
      if(fTGeoManager==nullptr){
        TGeoManager* old_global_manager = gGeoManager;
        gGeoManager = nullptr;
        fTGeoManager = new TGeoManager("NucDeEx","NucDeEx");
        fTGeoElementTable = fTGeoManager->GetElementTable();
        gGeoManager = old_global_manager;
        std::cout << "TGeoManager new" << std::endl;
      }
    }

    void SetPATH()
    {
      if(NUCDEEX_ROOT.length()>0) return;
      const char* env;
#ifndef WITH_NEUT
      env = std::getenv("NUCDEEX_ROOT");
      if(env==nullptr){
        std::cerr << "The env $NUCDEEX_ROOT is not specified" << std::endl;
        exit(1);
      }
      NUCDEEX_ROOT = (std::string) env;
#else
      env = std::getenv("NEUT_ROOT");
      if(env==nullptr){
        std::cerr << "The env $NEUT_ROOT is not specified" << std::endl;
        exit(1);
      }
      NUCDEEX_ROOT = (std::string) env + (std::string)"/share/nucdeex";
#endif
      std::cout << "NUCDEEX_ROOT: " << NUCDEEX_ROOT.c_str() << std::endl;
    }

#ifdef INCL_DEEXCITATION_NUCDEEX
    void SetPATH(G4INCL::Config* config)
    {
      if(NUCDEEX_ROOT.length()>0) return;
      NUCDEEX_ROOT = config->getNucDeExDataFilePath();
      std::cout << "NUCDEEX_ROOT: " << NUCDEEX_ROOT.c_str() << std::endl;
    }
#endif
  }
}
