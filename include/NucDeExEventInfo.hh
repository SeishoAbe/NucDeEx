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
  int EventID, fStatus, fShell;
  // fStatus:
  //     1: OK
  //     0: The nucleus is not supported
  //    -1: Error. Something wrong happens in energy conservation.
  int Zt, Nt, Z, N;
  double Ex;
  TVector3 Pinit; // Initial momentum of nucleus

  // particle level info
  std::vector<NucDeExParticle> ParticleVector;
};
#endif
