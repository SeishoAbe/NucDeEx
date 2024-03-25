#ifndef __NUCDEEXDEEXCITATIONTALYS__HH__
#define __NUCDEEXDEEXCITATIONTALYS__HH__

#include <string>
#include <ostream>
#include <sstream>
#include <vector>

#include "NucDeExConsts.hh"
#include "NucDeExDeexcitationBase.hh"

#include <TFile.h>
#include <TGraph.h>
#include <TParticle.h>

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

class NucDeExDeexcitationTALYS: public NucDeExDeexcitationBase{
  public:
  NucDeExDeexcitationTALYS();
  NucDeExDeexcitationTALYS(const int ld=1, const bool p_o=1);
#ifdef INCL_DEEXCITATION_NUCDEEX
  NucDeExDeexcitationTALYS(const int ld=1, const bool p_o=1, G4INCL::Config *config=0);
#endif
  virtual ~NucDeExDeexcitationTALYS(){};

  NucDeExEventInfo DoDeex(const int Zt, const int Nt,
                          const int Z, const int N, const double Ex,
                          const TVector3& mom=TVector3(0,0,0));
  // Zt, Nt  : Target nucleus Z and N (supports 16O or 12C currently)
  // Z, N    : Residual nucleus Z and N (having excitation energy)
  // Ex      : Excitation energy
  // mom     : 3D momentum of residual nucleus

  private:
  int DecayMode(const double Ex);

  // --- ROOT related methods & members --- //
  bool OpenROOT(const int Zt,const int Nt, const int Z, const int N);
  bool GetBrTGraph(const std::string st);
  int  GetBrExTGraph(const std::string st, const double ex_t, const int mode); 
  bool DaughterExPoint(double *d_Ex, int *d_point); //call by pointer
    // The nearest TGraph point will be returned
  void DeleteTGraphs();

  TFile* rootf;
  TGraph* g_br[NucDeEx::num_particle];
  TGraph* g_br_ex;

  int ldmodel;
  bool parity_optmodall;

};
#endif
