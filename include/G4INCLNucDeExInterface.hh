#ifndef G4INCLNucDeExInterface_hh
#define G4INCLNucDeExInterface_hh 1

#include "G4INCLConfig.hh"
#include "G4INCLIDeExcitation.hh"
#include "G4INCLEventInfo.hh"

#include "NucDeExDeexcitation.hh"

#include "TVector3.h"
#include <vector>

class G4INCLNucDeExInterface : public G4INCL::IDeExcitation {
public:
  G4INCLNucDeExInterface(G4INCL::Config*);
  virtual ~G4INCLNucDeExInterface();

  virtual void deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i);

private:
  G4INCL::Config *theConfig;
  NucDeExDeexcitation *theNucDeEx;
  int Zt,Nt,At;
  TVector3 Pinit;
  NucDeExEventInfo theNucDeExResult;
};

#endif
