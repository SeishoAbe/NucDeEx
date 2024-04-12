#include <iostream>
#include <string>
#include <cstdlib>

#include "NucDeExUtils.hh"

namespace NucDeEx{
  namespace Utils{
    // --- Define & initialize --- //
    int fVerbose=0;
    std::string NUCDEEX_ROOT="";
    NucDeExNucleusTable* NucleusTable=NULL;
    TDatabasePDG* fTDatabasePDG=NULL;
    TGeoManager* fTGeoManager=NULL;
    TGeoElementTable* fTGeoElementTable=NULL;

    void Init()
    {
      if(NucleusTable==NULL){
         NucleusTable = new NucDeExNucleusTable();
      }
      if(fTDatabasePDG==NULL){
        fTDatabasePDG = new TDatabasePDG();
      }
      if(fTGeoManager==NULL){
        fTGeoManager = new TGeoManager("NucDeEx","NucDeEx");
        fTGeoElementTable = fTGeoManager->GetElementTable();
      }
    }

    void SetPATH()
    {
      if(NUCDEEX_ROOT.length()>0) return;
      const char* env;
#ifndef WITH_NEUT
      env = std::getenv("NUCDEEX_ROOT");
      if(env==NULL){
        std::cerr << "The env $NUCDEEX_ROOT is not specified" << std::endl;
        exit(1);
      }
      NUCDEEX_ROOT = (std::string) env;
#else
      env = std::getenv("NEUT_ROOT");
      if(env==NULL){
        std::cerr << "The env $NEUT_ROOT is not specified" << std::endl;
        exit(1);
      }
      NUCDEEX_ROOT = (std::string) env + (std::string)"/../../src/nucdeex/nucdeex";
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
