#include <iostream>

#include "G4INCLNucDeExInterface.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLParticle.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLParticleSpecies.hh"
#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"

G4INCLNucDeExInterface::G4INCLNucDeExInterface(G4INCL::Config *config) :
  IDeExcitation(config),
  theConfig(config),
  theNucDeEx(new NucDeExDeexcitation(2,1,theConfig))
{
  NucDeEx::Utils::fVerbose=2;
  NucDeEx::Random::SetSeed(seed);
  Zt = theConfig->getTargetZ();
  At = theConfig->getTargetA();
  Nt = At-Zt;
}

G4INCLNucDeExInterface::~G4INCLNucDeExInterface() {
  delete theNucDeEx;
}

void G4INCLNucDeExInterface::deExciteRemnant(G4INCL::EventInfo *eventInfo, const int i) {
  int particleIndex = eventInfo->nParticles;

  Pinit.SetXYZ(eventInfo->pxRem[i],eventInfo->pyRem[i],eventInfo->pzRem[i]);

  theNucDeExResult = theNucDeEx->DoDeex(Zt,Nt,
                     eventInfo->ZRem[i],eventInfo->ARem[i]-eventInfo->ZRem[i],
                     0,eventInfo->EStarRem[i],Pinit);
  std::vector<NucDeExParticle> ParticleVector = theNucDeExResult.ParticleVector;
  int size=ParticleVector.size();

  for(int j = 0; j < size ; ++j) { // Copy NucDeEx result to the EventInfo
    NucDeExParticle particle = ParticleVector.at(j);
    if(!particle._flag) continue; // skip intermediate states
    const int PDG = particle._PDG;
    if(PDG==22){// gamma
      eventInfo->A[particleIndex] = 0;
      eventInfo->Z[particleIndex] = 0;
      eventInfo->PDGCode[particleIndex] = PDG;
    }else if(PDG==2112){ // neutron
      eventInfo->A[particleIndex] = 1;
      eventInfo->Z[particleIndex] = 0;
      eventInfo->PDGCode[particleIndex] = PDG;
    }else if(PDG==2212){ // proton
      eventInfo->A[particleIndex] = 1;
      eventInfo->Z[particleIndex] = 1;
      eventInfo->PDGCode[particleIndex] = PDG;
    }else{
      eventInfo->A[particleIndex] = (PDG%10000)/10;
      eventInfo->Z[particleIndex] = (PDG%10000000-eventInfo->A[particleIndex]*10)/10000;
      eventInfo->PDGCode[particleIndex] = eventInfo->Z[particleIndex]*1000 + eventInfo->A[particleIndex];
    }
    //std::cout << "PDG = " << PDG << "   Z = " << eventInfo->Z[particleIndex]
    //          << "   A = " << eventInfo->A[particleIndex] << std::endl;

    eventInfo->S[particleIndex] = 0; // not supported
    eventInfo->emissionTime[particleIndex] = -1.0;
    eventInfo->EKin[particleIndex] = particle.kE();
    const double px = particle._momentum.X();
    const double py = particle._momentum.Y();
    const double pz = particle._momentum.Z();
    const double plab = particle._momentum.Mag();
    eventInfo->px[particleIndex] = px;
    eventInfo->py[particleIndex] = py;
    eventInfo->pz[particleIndex] = pz;
    double pznorm =0.;
    if(plab>0.)pznorm = pz/plab;
    eventInfo->theta[particleIndex] = G4INCL::Math::toDegrees(G4INCL::Math::arcCos(pznorm));
    eventInfo->phi[particleIndex] = G4INCL::Math::toDegrees(std::atan2(py,px));
    eventInfo->origin[particleIndex] = i+0; // Origin: De-excitation (not supported)
    eventInfo->history.push_back("cxx"); // No history tracking for NucDeEx
    eventInfo->ParticleBias[particleIndex] = G4INCL::Particle::getTotalBias();
    //
    particleIndex++;
  }

  eventInfo->nParticles = particleIndex;
}

