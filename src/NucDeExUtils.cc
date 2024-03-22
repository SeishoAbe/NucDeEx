#include <iostream>
#include <string>

#include "NucDeExUtils.hh"

int NucDeExUtils::fVerbose=0;
int NucDeExUtils::fSeed=0;
std::string NucDeExUtils::NUCDEEX_ROOT="";

/////////////////////////////////////////////
void NucDeExUtils::SetPATH()
/////////////////////////////////////////////
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
/////////////////////////////////////////////
void NucDeExUtils::SetPATH(G4INCL::Config* config)
/////////////////////////////////////////////
{
  if(NUCDEEX_ROOT.length()>0) return;
  NUCDEEX_ROOT = config->getNucDeExDataFilePath();
  std::cout << "NUCDEEX_ROOT: " << NUCDEEX_ROOT.c_str() << std::endl;
}
#endif
