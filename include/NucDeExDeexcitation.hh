#ifndef __NUCDEEXDEEXCITATION__HH__
#define __NUCDEEXDEEXCITATION__HH__

#include <string>
#include <ostream>
#include <sstream>
#include <vector>

#include "NucDeExConsts.hh"
#include "NucDeExDeexcitationBase.hh"
#include "NucDeExDeexcitationTALYS.hh"
#include "NucDeExDeexcitationPhole.hh"

#include <TFile.h>
#include <TGraph.h>
#include <TParticle.h>

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

class NucDeExDeexcitation: public NucDeExDeexcitationBase{
  public:
  NucDeExDeexcitation();
  NucDeExDeexcitation(const int ld=1, const bool p_o=1);
#ifdef INCL_DEEXCITATION_NUCDEEX
  NucDeExDeexcitation(const int ld=1, const bool p_o=1, G4INCL::Config *config=0);
#endif
  // ld: ldmodel, p_o: parity_optmodall

  virtual ~NucDeExDeexcitation();

  // --- Main function ---//
  NucDeExEventInfo DoDeex(const int Zt, const int Nt,
                          const int Z, const int N, const int shell, const double Ex,
                          const TVector3& mom=TVector3(0,0,0));
  // Zt, Nt  : Target nucleus Z and N (supports 16O or 12C currently)
  // Z, N    : Residual nucleus Z and N (having excitation energy)
  // shell  0 : Shell level will be determined according to Ex (box cut)
  //        1 : s1/2-hole
  //        2 : p3/2-hole
  //        3 : p1/2-hole (only for 16O target)
  // Ex      : Only used if shell==1
  // mom     : 3D momentum of residual nucleus
  // 
  // return  0 : The target or residual nucleus is not supported
  //           : No root file 
  //           : The g.s. nucleus is added in this case.
  //         1 : Success
  //        -1 : Fatal error 
  //           : Strange target & residual nuclei not in nucleus table, no mass profile
  //           : no particle info is added
  
  int ExtoShell(const int Zt, const int Nt, const double Ex);
  int GetShell(){return fShell;};

  private:
  NucDeExDeexcitationTALYS* deex_talys;
  NucDeExDeexcitationPhole* deex_phole;

  // --- model parameters --- //
  int ldmodel;
  bool parity_optmodall;

  // --- for Ex to shell --- //
  int fShell;
  const double Ex_12C_s12=16.0;
  const double Ex_16O_s12=16.0;
  const double Ex_16O_p32=4.0;
};
#endif
