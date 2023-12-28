//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
#define NUCDEEX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4NucDeExInterface_hh
#define G4NucDeExInterface_hh 1

#ifdef NUCDEEX_IN_GEANT4_MODE

#include "G4VPreCompoundModel.hh"
#include "G4ReactionProduct.hh"
#include "G4Fragment.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

#include "Deexcitation.hh"

#include "TVector3.h"
#include <vector>

class G4NucDeExInterface : public G4VPreCompoundModel {
public:
  G4NucDeExInterface();
  virtual ~G4NucDeExInterface();

  virtual G4ReactionProductVector *DeExcite(G4Fragment &aFragment);

  virtual G4HadFinalState *ApplyYourself(G4HadProjectile const &, G4Nucleus &) {
    return NULL;
  }

  virtual void ModelDescription(std::ostream& outFile) const;
  virtual void DeExciteModelDescription(std::ostream& outFile) const;

private:
  Deexcitation *theNucDeEx;
  int Zt,Nt,At;
  TVector3 Pinit;

  std::vector<Particle> *theNucDeExResult;
  G4long eventNumber;

  /// \brief Convert an NucDeEx particle to a G4DynamicParticle
  G4ReactionProduct *toG4Particle(G4int A, G4int Z, G4int S, G4double kinE, G4double px, G4double py, G4double pz) const;

  /// \brief Convert A, Z and S to a G4ParticleDefinition
  G4ParticleDefinition *toG4ParticleDefinition (G4int A, G4int Z, G4int S) const;

};

#endif // NUCDEEX_IN_GEANT4_MODE

#endif
