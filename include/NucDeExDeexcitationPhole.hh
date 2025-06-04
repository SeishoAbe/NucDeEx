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
  NucDeExDeexcitationPhole(int v);
    // 11B* and 11C* are common (Panin et al., Phys. Lett. B 753 204-210. Experimental data)
    // 15N* and 15O* are as follows:
    //  0 -> Based on Ejili, Phys. Rev. C 58, 3
    //  1 -> Based on Yosoi (D-thesis)
  ~NucDeExDeexcitationPhole(){};

  void SetParameters();

  NucDeExEventInfo DoDeex(const int Zt, const int Nt,
                          const int Z, const int N, const double Ex,
                          const TVector3& mom=TVector3(0,0,0));
    // Note: Energy level is determined irrelevant to Ex


  void DoDeex_v3(const int Zt, const int Nt,
                 const int Z, const int N, const double Ex,
                 const TVector3& mom=TVector3(0,0,0));
  // Called if version==3
  //   - This function is expected to be used with new carbon SF, PRC 110, 054612 (2024).
  //   - We don't need to determine the excited state with this new SF.

  void DoDeex_v2(const int Zt, const int Nt,
                 const int Z, const int N, const double Ex,
                 const TVector3& mom=TVector3(0,0,0));
  // Called if version==2
  //   - p1/2 g.s. is included, so it should be distinguished in this function
  //   - The excited states above the separation energy is already excluded in the upstream process.
  //     They are handled by NucDeExDeexcitationTALYS.

  void DoDeex_v1(const int Zt, const int Nt,
                 const int Z, const int N, const double Ex,
                 const TVector3& mom=TVector3(0,0,0));
  // Called if version==1
  //   - p1/2 g.s. is already excluded in upstream process (NucDeExDeexcitation)
  //   - The excited states above the separation energy is included
  //   - The excited state are decided **irrelevant** to the inputed excitation energy

  void SetPtrTALYS(NucDeExDeexcitationTALYS* t){deex_talys = t;};

  private:
  int version;
  NucDeExDeexcitationTALYS* deex_talys;

  int Nlevel_p32_11B;
  double *E_p32_11B, *Br_p32_11B;

  int Nlevel_p32_11C;
  double *E_p32_11C, *Br_p32_11C;

  int Nlevel_p32_15N;
  double *E_p32_15N, *Br_p32_15N;

  int Nlevel_p32_15O;
  double *E_p32_15O, *Br_p32_15O;
};
#endif
