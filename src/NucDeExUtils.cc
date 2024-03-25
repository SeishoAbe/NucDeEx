#include <iostream>
#include <string>

#include "NucDeExUtils.hh"

namespace NucDeEx{
  namespace Utils{
    int fVerbose=0;
    std::string NUCDEEX_ROOT="";
    NucDeExNucleusTable* NucleusTable=NULL;

    void SetPATH()
    {
      if(NUCDEEX_ROOT.length()>0) return;
      const char* env = getenv("NUCDEEX_ROOT");
      NUCDEEX_ROOT = env;
#ifdef WITH_NEUT
      env = getenv("NEUT_ROOT");
      NUCDEEX_ROOT = env + (std::string)"/../../src/nucdeex/nucdeex";
#endif
      if(env==NULL){
        std::cerr << "PATH to nucleus table is not specified" << std::endl;
        exit(1);
      }
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
