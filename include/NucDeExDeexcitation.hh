#ifndef __NUCDEEXDEEXCITATION__HH__
#define __NUCDEEXDEEXCITATION__HH__

#include <string>
#include <ostream>
#include <sstream>
#include <vector>

#include "NucDeExDeexcitationBase.hh"
#include "NucDeExDeexcitationTALYS.hh"
#include "NucDeExDeexcitationPhole.hh"

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

class NucDeExDeexcitation: public NucDeExDeexcitationBase{
  public:
  NucDeExDeexcitation();
  NucDeExDeexcitation(const int ld, const bool p_o, const int v);
#ifdef INCL_DEEXCITATION_NUCDEEX
  NucDeExDeexcitation(const int ld, const bool p_o, const int v, G4INCL::Config *config);
#endif
  // ld: ldmodel, p_o: parity_optmodall

  virtual ~NucDeExDeexcitation();

  void Init();

  // --- Main function ---//
  NucDeExEventInfo DoDeex(const int Zt, const int Nt, const int Z, const int N,
                          const double Ex, const TVector3& mom=TVector3(0,0,0));
  // New main method from v2.1. This should be used with version_phole=2
  // Zt, Nt  : Target nucleus Z and N (supports 16O or 12C currently)
  // Z, N    : Residual nucleus Z and N (having excitation energy)
  // Ex      : Excitation energy
  // mom     : 3D momentum of residual nucleus

  NucDeExEventInfo DoDeex(const int Zt, const int Nt, const int Z, const int N,
                          const int shell, const double Ex,
                          const TVector3& mom=TVector3(0,0,0));
  // Old main method until v1.3. This should be used with version_phole=1
  // shell  0 : Shell level will be determined according to Ex (box cut)
  //        1 : s1/2-hole
  //        2 : p3/2-hole
  //        3 : p1/2-hole (only for 16O target)

  void SetVersionPhole(int v){ version_phole = v ; }
  int GetVersionPhole(){ return version_phole ; }

  private:
  NucDeExDeexcitationTALYS* deex_talys;
  NucDeExDeexcitationPhole* deex_phole;

  // --- model parameters --- //
  int ldmodel;
  bool parity_optmodall;
  int version_phole;

  // --- for Ex to shell (for old methond until v1.3) --- //
  int ExtoShell(const int Zt, const int Nt, const double Ex);
  int GetShell(){return fShell;};
  int fShell;
  const double Ex_12C_s12=16.0;
  const double Ex_16O_s12=16.0;
  const double Ex_16O_p32=4.0;
};
#endif
