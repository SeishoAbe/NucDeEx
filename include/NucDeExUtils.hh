#ifndef __NUCDEEXUTILS__HH__
#define __NUCDEEXUTILS__HH__

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif


class NucDeExUtils{
  public:
  NucDeExUtils(){;}
  ~NucDeExUtils(){;}

  static void SetVerbose(int v){ fVerbose = v; };
  static int  GetVerbose(){ return fVerbose; };
  static void SetSeed(int v){ fSeed = v; };
  static int  GetSeed(){ return fSeed; };

  static void SetPATH();
#ifdef INCL_DEEXCITATION_NUCDEEX
  static void SetPATH(G4INCL::Config* config);
#endif
  static std::string GetPATH(){ return NUCDEEX_ROOT; };

  private:
  static int fVerbose;
  static int fSeed;
  static std::string NUCDEEX_ROOT;
  // Equivalent to env $NUCDEEX_ROOT for standalone use
};
#endif
