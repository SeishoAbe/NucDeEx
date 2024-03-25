#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "NucDeExParticle.hh"
#include "NucDeExDeexcitationBase.hh"

///////////////////////////
NucDeExDeexcitationBase::NucDeExDeexcitationBase(): EventID(0)
///////////////////////////
{
  if(NucDeEx::Utils::NucleusTable==NULL){
    NucDeEx::Utils::NucleusTable = new NucDeExNucleusTable();
  }
  fTDatabasePDG = new TDatabasePDG();
  fTGeoManager = new TGeoManager("NucDeEx","NucDeEx");
  fTGeoElementTable = fTGeoManager->GetElementTable();
}

///////////////////////////
NucDeExDeexcitationBase::~NucDeExDeexcitationBase()
///////////////////////////
{
  delete fTDatabasePDG;
  delete fTGeoManager;
  //delete fTGeoElementTable;
}

///////////////////////////
void NucDeExDeexcitationBase::Init()
///////////////////////////
{
  NucDeEx::Utils::NucleusTable->ReadTables(0);
}

/////////////////////////////////////////////
void NucDeExDeexcitationBase::SaveEventLevelInfo(const int Zt, const int Nt,
                        const int Z, const int N,const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
  // event level info
  EventInfo.EventID = EventID;
  EventInfo.Zt     = Zt;
  EventInfo.Nt     = Nt;
  EventInfo.Z      = Z;
  EventInfo.N      = N;
  EventInfo.Ex     = Ex;
  EventInfo.Pinit  = mom;
}

/////////////////////////////////////////////
bool NucDeExDeexcitationBase::Decay(const bool breakflag)
/////////////////////////////////////////////
{
  if(NucDeEx::Utils::fVerbose>0){
    std::cout << "NucDeExDeexcitationBase::Decay()" << std::endl;
  }
  bool status=1; // sucess

  // ---- CM frame --- //
  // --- Calculate kinematics (CM)
  //    mass used in the following calculation should include excitation energy
  double mass_ex_daughter = mass_daughter + Ex_daughter;

  double cmMomentum = std::sqrt(Qvalue*(Qvalue + 2.*mass_particle)*
                              (Qvalue + 2.*mass_ex_daughter)*
                              (Qvalue + 2.*mass_particle + 2.*mass_ex_daughter) )/
                              (Qvalue + mass_particle + mass_ex_daughter)/2.;
  double kE_particle = sqrt( pow(cmMomentum,2) + pow(mass_particle,2) ) - mass_particle;
  double kE_daughter = sqrt( pow(cmMomentum,2) + pow(mass_ex_daughter,2) ) - mass_ex_daughter;
  double kE_sum = kE_particle + kE_daughter; // just for check
  if( Qvalue>0 && (kE_sum - Qvalue)/Qvalue > check_criteria){
    std::cerr << "WARNING @ NucDeExDeexcitationBase::Decay: Something wrong in kinematics calculation" << std::endl;
    status=0;
  }
  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "mass_target   = " << mass_target << std::endl;
    std::cout << "mass_particle = " << mass_particle << std::endl;
    std::cout << "mass_daughter = " << mass_daughter << std::endl;
    std::cout << "cmMomentum  = " << cmMomentum << std::endl;
    std::cout << "kE_particle = " << kE_particle << std::endl;
    std::cout << "kE_daughter = " << kE_daughter << std::endl;
    std::cout << "kE_sum      = " << kE_sum << " <- consistent with Qvalue = " << Qvalue << std::endl;
  }
  
  // --- Detemine momentum direction (uniform) (CM)
  // --- Then set momentum vectors
  double costheta = 2.*NucDeEx::Random::random()-1.; // [-1, 1]
  double sintheta = sqrt( 1. - pow(costheta,2) );
  double phi      = 2*TMath::Pi()*NucDeEx::Random::random(); // [0,2pi]
  TVector3 dir( sintheta*cos(phi), sintheta*sin(phi), costheta );
  mom_particle = -1*cmMomentum*dir; // -P
  mom_daughter = cmMomentum*dir; // P
  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "dir:          "; dir.Print();
    std::cout << "mom_particle: ";mom_particle.Print();
    std::cout << "mom_daughter: ";mom_daughter.Print();
  }
  

  // ---- CM -> LAB --- //
  // --- Get total energy of target (i.e., CM frame) in the LAB frame for boost
  //     This is "Total CM energy" (in the LAB frame)
  //     This should not include excitation energy
  double totalE_target = sqrt( pow(mom_target.Mag(),2) + pow(mass_target,2) );

  // --- Store info as Particle. then boost it
  NucDeExParticle p_particle(NucDeEx::PDG_particle[decay_mode],
                      mass_particle,
                      mom_particle,
                      NucDeEx::particle_name[decay_mode],
                      1,0);// trace flag==1, Ex_daughter==0
  double totalE_particle_bef = p_particle.totalE();
  p_particle.Boost(totalE_target,mom_target);// BOOST!
  double totalE_particle_aft = p_particle.totalE();

  NucDeExParticle p_daughter(PDGion(Z_daughter,N_daughter),
                      mass_daughter, // w/o excitation E
                      mom_daughter,
                      name_daughter,
                      0,Ex_daughter); // intermediate state in default

  double totalE_daughter_bef = p_daughter.totalE();
  p_daughter.Boost(totalE_target,mom_target);
  double totalE_daughter_aft = p_daughter.totalE();

  // 
  double kE_target = totalE_target - mass_target;
  double totalE_ex_target = totalE_target + Ex_target; // w/ excitation E

  double totalE_bef = totalE_particle_bef + totalE_daughter_bef;
  double totalE_aft = totalE_particle_aft + totalE_daughter_aft;
  double totalE_ex_bef = totalE_bef + Ex_daughter; // w/ excitation E
  double totalE_ex_aft = totalE_aft + Ex_daughter;

  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "totalE_target = " << totalE_target << std::endl;
    std::cout << "kE_target = " << kE_target << std::endl;
    std::cout << "Ex_target = " << Ex_target << std::endl;
    std::cout << "totalE_ex_target = " << totalE_ex_target << std::endl;
    std::cout << "totalE_bef = " << totalE_bef << std::endl;
    std::cout << "totalE_aft = " << totalE_aft << std::endl;
    std::cout << "  diff = " << totalE_aft-totalE_bef << std::endl;
    std::cout << "totalE_ex_bef = " << totalE_ex_bef << std::endl;
    std::cout << "totalE_ex_aft = " << totalE_ex_aft << std::endl;
  }


  // --- Check Energy conservation (4dim energy including momentum)
  //       Fundamental energy conservation 
  //       (Total energy in LAB) = (Total energy in CM after boost)
  // i.e., (Total energy of target w/ ex in LAB) 
  //         = (Total energy of particle after boost) + (Total energy of daughter w/ ex after boost)
  //            The last two terms are calculated from CM at first, and then boosted.
  if(totalE_ex_target>0 && (totalE_ex_aft-totalE_ex_target)/totalE_ex_target>check_criteria){
    std::cerr << "WARNING: @ NucDeExDeexcitation:Decay(): Energy is not conserved..." << std::endl;
    status=0;
  }

  // Then, push back
  EventInfo.ParticleVector.push_back(p_particle);
  // DoDeex loop will be end -> turn on track flag, because it is not intermediate state
  if(breakflag) p_daughter._flag=1;
  EventInfo.ParticleVector.push_back(p_daughter);

  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "flag_particle = " << p_particle._flag << std::endl;
    std::cout << "flag_daughter = " << p_daughter._flag << std::endl;
  }
  return status;
}

/////////////////////////////////////////////
void NucDeExDeexcitationBase::AddGSNucleus(const int Z,const int N, const TVector3& mom)
/////////////////////////////////////////////
{
  double mass_target = ElementMassInMeV(Z+N, Z);
  NucDeExNucleus* nuc_target = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z,N);
  if(nuc_target==NULL || mass_target<0){
    EventInfo.fStatus=0;
    return; // no particle is added to the vector
  }
  string name_target = (string)nuc_target->name;
  NucDeExParticle nucleus(PDGion(Z,N),
                          mass_target, // w/o excitation E
                          mom,
                          name_target, 
                          1,0); // track flag on // zero ex
  if(NucDeEx::Utils::fVerbose>0){
    std::cout << "AddGSNucleus(): " << name_target << std::endl;
  }
  EventInfo.fStatus=1;
  EventInfo.ParticleVector.push_back(nucleus);
  return;
}

/////////////////////////////////////////////
const double NucDeExDeexcitationBase::ElementMassInMeV(const int A, const int Z)
/////////////////////////////////////////////
{
  return ElementMassInMeV( fTGeoElementTable->GetElementRN(A,Z) );
}

/////////////////////////////////////////////
const double NucDeExDeexcitationBase::ElementMassInMeV(const TGeoElementRN* ele)
/////////////////////////////////////////////
{
  if(ele==0) return -1; // no profile can be found.
  double mass = ele->MassNo()*NucDeEx::amu_c2
                   + ele->MassEx(); // (MeV)
  double mass_amu = mass/NucDeEx::amu_c2;
  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "mass (MeV) = " << mass 
         << "    (amu) = " << mass_amu << std::endl;
  }
  return mass;
}

/////////////////////////////////////////////
const int NucDeExDeexcitationBase::PDGion(const int Z, const int N)
/////////////////////////////////////////////
{
  int pdg_int= 1e9 + Z*1e4 + (Z+N)*1e1;
  return pdg_int;
}
