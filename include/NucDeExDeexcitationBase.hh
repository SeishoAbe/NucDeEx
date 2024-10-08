#ifndef __NUCDEEXDEEXCITATIONBASE__HH__
#define __NUCDEEXDEEXCITATIONBASE__HH__

#include <string>
#include <ostream>
#include <sstream>
#include <vector>

#include "NucDeExEventInfo.hh"
#include "NucDeExNucleus.hh"
#include <TGeoElement.h>

#ifdef INCL_DEEXCITATIONBASE_NUCDEEX
#include "G4INCLConfig.hh"
#endif

class NucDeExDeexcitationBase{
  public:
  NucDeExDeexcitationBase();
  virtual ~NucDeExDeexcitationBase();

  protected:
  void SaveEventLevelInfo(const int Zt, const int Nt,
                          const int Z, const int N,const double Ex, 
                          const TVector3& mom);

  void Decay(const bool breakflag);
    // Calculate two-body decay considering Lorentz boost
    // If something wrong happens in energy conservation, set `EventInfo.fStatus=-1`.
  
  void AddGSNucleus(const int Z, const int N, const TVector3& mom=TVector3(0,0,0));
    // Add g.s. nucleus to EventInfo
    // If it is not listed in nucleustable, set `EventInfo.fStatus=0`.
  double ElementMassInMeV(const int A, const int Z);
  double ElementMassInMeV(const TGeoElementRN* ele);
  int PDGion(const int Z, const int N); // Z, N to PDG

  //--- params --- //
  // target nucleus info
  int Z_target, N_target;
  double Ex_target;
  double mass_target;
  TVector3 mom_target;
  NucDeExNucleus* nuc_target;
  std::string name_target;
  // daughter nucleus info
  int Z_daughter, N_daughter;
  double Ex_daughter;
  double mass_daughter;
  TVector3 mom_daughter;
  NucDeExNucleus* nuc_daughter;
  std::string name_daughter;
  // decay particle info
  double mass_particle;
  TVector3 mom_particle;
  // decay info
  int decay_mode;
  double S;
  double Qvalue;

  // --- scoring --- //
  NucDeExEventInfo EventInfo;
  int EventID;
  
  // --- utils --- //
  std::ostringstream os;
  const double check_criteria=5e-3;
};
#endif
