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
  // fShell:
  //  v2.1~
  //    0: Did nothing
  //    1: TALYS (NucDeExDeexcitationTALYS)
  //    2: gamma discrete
  //    3: g.s. (AddGSNucleus)
  //  ~v1.3
  //    1 : s1/2-hole
  //    2 : p3/2-hole
  //    3 : p1/2-hole (only for 16O target)
  int Zt, Nt, Z, N;
  double Ex;
  TVector3 Pinit; // Initial momentum of nucleus

  // particle level info
  std::vector<NucDeExParticle> ParticleVector;
};
#endif
