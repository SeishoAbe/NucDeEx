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

#ifdef NUCDEEX_IN_GEANT4_MODE

#include "G4NucDeExInterface.hh"
#include "G4ParticleDefinition.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <cmath>

G4NucDeExInterface::G4NucDeExInterface() :
  G4VPreCompoundModel(NULL, "NucDeEx"),
  theNucDeEx(new NucDeExDeexcitation(2,1)),
  eventNumber(0)
{
  theNucDeEx->SetSeed(1);
  theNucDeEx->SetVerbose(0);
}

G4NucDeExInterface::~G4NucDeExInterface() {
  delete theNucDeEx;
}

G4ReactionProductVector *G4NucDeExInterface::DeExcite(G4Fragment &aFragment) {
  const G4int ARem = aFragment.GetA_asInt();
  const G4int ZRem = aFragment.GetZ_asInt();
  const G4double eStarRem = aFragment.GetExcitationEnergy() / MeV;
  //const G4double jRem = aFragment.GetAngularMomentum().mag() / hbar_Planck; // unused
  const G4LorentzVector &pRem = aFragment.GetMomentum();
  const G4double pxRem = pRem.x() / MeV;
  const G4double pyRem = pRem.y() / MeV;
  const G4double pzRem = pRem.z() / MeV;
  Pinit.SetXYZ(pxRem,pyRem,pzRem);

  eventNumber++;

  theNucDeEx->DoDeex(Zt,Nt,
                     ZRem,ARem-ZRem, 
                     0,eStarRem,Pinit);

  G4ReactionProductVector *result = new G4ReactionProductVector;

  theNucDeExResult = theNucDeEx->GetParticleVector();
  int size = theNucDeExResult->size();

  for(int j = 0; j < size ; ++j) { // Copy NucDeEx result to the EventInfo
    NucDeExParticle particle = theNucDeExResult->at(j);
    if(!particle._flag) continue; // skip intermediate states
    const int PDG = particle._PDG;
    int A,Z;
    if(PDG==22){// gamma
      A = 0;
      Z = 0;
    }else if(PDG==2112){ // neutron
      A = 1;
      Z = 0;
    }else if(PDG==2212){ // proton
      A = 1;
      Z = 1;
    }else{
      A = (PDG%10000)/10;
      Z = (PDG%10000000-A*10)/10000;
    }

    G4ReactionProduct *product = toG4Particle(A,Z,0, // S
                                              particle.kE(),
                                              particle._momentum.X(),
                                              particle._momentum.Y(),
                                              particle._momentum.Z());

    if(product)
      result->push_back(product);
  }
  return result;
}

G4ParticleDefinition *G4NucDeExInterface::toG4ParticleDefinition(G4int A, G4int Z, G4int S) const {
  if     (A == 1 && Z == 1 && S == 0)  return G4Proton::Proton();
  else if(A == 1 && Z == 0 && S == 0)  return G4Neutron::Neutron();
  else if(A == 1 && Z == 0 && S == -1)  return G4Lambda::Lambda();
  else if(A == -1 && Z == 1 && S == 0)  return G4PionPlus::PionPlus();
  else if(A == -1 && Z == -1 && S == 0) return G4PionMinus::PionMinus();
  else if(A == -1 && Z == 0 && S == 0)  return G4PionZero::PionZero();
  else if(A == 0 && Z == 0 && S == 0)  return G4Gamma::Gamma();
  else if(A == 2 && Z == 1 && S == 0)  return G4Deuteron::Deuteron();
  else if(A == 3 && Z == 1 && S == 0)  return G4Triton::Triton();
  else if(A == 3 && Z == 2 && S == 0)  return G4He3::He3();
  else if(A == 4 && Z == 2 && S == 0)  return G4Alpha::Alpha();
  else if(A > 0 && Z > 0 && A > Z) { // Returns ground state ion definition.
    return G4IonTable::GetIonTable()->GetIon(Z, A, std::abs(S));//S is the number of lambdas
  } else { // Error, unrecognized particle
    G4cout << "Can't convert particle with A=" << A << ", Z=" << Z << ", S=" << S << " to G4ParticleDefinition, trouble ahead" << G4endl;
    return 0;
  }
}

G4ReactionProduct *G4NucDeExInterface::toG4Particle(G4int A, G4int Z, G4int S,
						 G4double kinE,
						 G4double px,
                                                 G4double py, G4double pz) const {
  G4ParticleDefinition *def = toG4ParticleDefinition(A, Z, S);
  if(def == 0) { // Check if we have a valid particle definition
    return 0;
  }

  const G4double energy = kinE * MeV;
  const G4ThreeVector momentum(px, py, pz);
  const G4ThreeVector momentumDirection = momentum.unit();
  G4DynamicParticle p(def, momentumDirection, energy);
  G4ReactionProduct *r = new G4ReactionProduct(def);
  (*r) = p;
  return r;
}

void G4NucDeExInterface::ModelDescription(std::ostream& outFile) const {
   outFile << "NUCDEEX++ does not provide an implementation of the ApplyYourself method!\n\n";
}

void G4NucDeExInterface::DeExciteModelDescription(std::ostream& outFile) const {
   outFile 
     << "NUCDEEX++ is a statistical model for nuclear de-excitation. It simulates\n"
     << "the gamma emission and the evaporation of neutrons, light charged particles\n"
     << "and IMFs, as well as fission where applicable. The code included in Geant4\n"
     << "is a C++ translation of the original Fortran code NUCDEEX07. Although the model\n"
     << "has been recently extended to hypernuclei by including the evaporation of lambda\n"
     << "particles. More details about the physics are available in the\n"
     << "Geant4 Physics Reference Manual and in the reference articles.\n\n"
     << "References:\n"
     << "(1) A. Kelic, M. V. Ricciardi, and K. H. Schmidt, in Proceedings of Joint\n"
     << "ICTP-IAEA Advanced Workshop on Model Codes for Spallation Reactions,\n"
     << "ICTP Trieste, Italy, 4–8 February 2008, edited by D. Filges, S. Leray, Y. Yariv,\n"
     << "A. Mengoni, A. Stanculescu, and G. Mank (IAEA INDC(NDS)-530, Vienna, 2008), pp. 181–221.\n\n"
     << "(2) J.L. Rodriguez-Sanchez, J.-C. David et al., Phys. Rev. C 98, 021602 (2018)\n\n"; 
}

#endif // NUCDEEX_IN_GEANT4_MODE
