#ifndef __NUCDEEXEVENTINFO__HH__
#define __NUCDEEXEVENTINFO__HH__

#include <vector>

#include "NucDeExParticle.hh"

class NucDeExEventInfo{
  public:
  NucDeExEventInfo();
  ~NucDeExEventInfo();

  void InitParameters();

  // event level info
  int eventID, fStatus, fShell;
  int Zt, Nt, Z, N;
  double Ex, S, MissE;
  TVector3 Pinit; // Initial momentum of nucleus
  // particle level info
  std::vector<NucDeExParticle> ParticleVector;
};
#endif
