#ifndef __NUCDEEXDEEXCITATIONPHOLE__HH__
#define __NUCDEEXDEEXCITATIONPHOLE__HH__

#include <string>
#include <ostream>
#include <sstream>
#include <vector>

#include "NucDeExConsts.hh"
#include "NucDeExDeexcitationBase.hh"
#include "NucDeExDeexcitationTALYS.hh"

class NucDeExDeexcitationPhole: public NucDeExDeexcitationBase{
  public:
  NucDeExDeexcitationPhole();
  NucDeExDeexcitationPhole(int f);
    // 11B* and 11C* are common (Panin et al., Phys. Lett. B 753 204-210. Experimental data)
    // 15N* and 15O* are as follows:
    //  0 -> Based on Ejili, Phys. Rev. C 58, 3
    //  1 -> Based on Yosoi (D-thesis)
  ~NucDeExDeexcitationPhole(){};

  NucDeExEventInfo DoDeex(const int Zt, const int Nt,
                          const int Z, const int N, const double Ex,
                          const TVector3& mom=TVector3(0,0,0));
    // Note: Energy level is determined irrelevant to Ex


  void SetPtrTALYS(NucDeExDeexcitationTALYS* t){deex_talys = t;};

  private:
  int flag_model;
  NucDeExDeexcitationTALYS* deex_talys;
  

  // Constants for (p3/2)-1 Br
  // 11B* (Panin et al., Phys. Lett. B 753 204-210. Experimental data)
  static const int Nlevel_p32_11B = 3; 
  const double E_p32_11B[Nlevel_p32_11B]={0.,2.125,5.020};
  const double Br_p32_11B[Nlevel_p32_11B]={0.82,0.10,0.08};
  // 11C* (assume analogy of 11B*)
  //      the same Br, but energy is different
  static const int Nlevel_p32_11C = 3;
  const double E_p32_11C[Nlevel_p32_11C]={0.,2.000,4.804};
  const double Br_p32_11C[Nlevel_p32_11C]={0.82,0.10,0.08};

  // 15N* (Ejili, Phys. Rev. C 58, 3)
  //     Deexcitation from 9.93 MeV will be described from TALYS data
  //     Deexcitation from 10.70 is taken from Ejiri (100% proton)
  static const int Nlevel_p32_15N = 3;
  const double E_p32_15N[Nlevel_p32_15N]={6.32,9.93,10.70};
  const double Br_p32_15N[Nlevel_p32_15N]={0.872,0.064,0.064};
  // 15O* (Ejiri , Phys. Rev. C 58, 3)
  static const int Nlevel_p32_15O=3;
  const double E_p32_15O[Nlevel_p32_15O]={6.18,9.61,10.48};
  const double Br_p32_15O[Nlevel_p32_15O]={0.872,0.064,0.064}; // guess
  //const double Br_p32_15O[Nlevel_p32_15O]={1.,0,0}; // original Ejiri's value
};
#endif
