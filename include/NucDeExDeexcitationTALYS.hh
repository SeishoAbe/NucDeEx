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
  NucDeExDeexcitationTALYS(const int ld, const bool p_o);
#ifdef INCL_DEEXCITATION_NUCDEEX
  NucDeExDeexcitationTALYS(const int ld, const bool p_o, G4INCL::Config *config);
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
  void GetAllTGraph();
  std::map<std::string, int> map_root; 
    // store nucleus of root file: nucleus name -> [file]
  std::map<std::string, int> map_nucleus[NucDeEx::bins]; 
    // store nucleus of tgraph (target nucleus): nucleus name -> [nucleus]
  std::map<std::string, int> :: iterator it;
  int getRootID(const char* name);
  int getNucleusID(int i, const char* name);
  int fRootID, fNucleusID, fPoint;
    // [file], [nucleus], [Exbin]
    // use "decay_mdoe for [particle]

  int  GetBrExTGraph(const std::string st, const double ex_t, const int mode); 
  bool DaughterExPoint(double *d_Ex, int *d_point); //call by pointer

  TGraph* g_br_all[NucDeEx::bins][NucDeEx::bins][NucDeEx::num_particle]; 
    // [file][nucleus][particle]
  TGraph* g_br_ex_all[NucDeEx::bins][NucDeEx::bins][NucDeEx::num_particle][NucDeEx::bins];
    // [file][nucleus][particle][Exbin]


  int ldmodel;
  bool parity_optmodall;

};
#endif
