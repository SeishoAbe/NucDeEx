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
#include "NucDeExDeexcitationPhole.hh"


///////////////////////////
NucDeExDeexcitationPhole::NucDeExDeexcitationPhole(): flag_model(0)
///////////////////////////
{}
///////////////////////////
NucDeExDeexcitationPhole::NucDeExDeexcitationPhole(int f): flag_model(f)
///////////////////////////
{}


/////////////////////////////////////////////
NucDeExEventInfo NucDeExDeexcitationPhole::DoDeex(const int Zt, const int Nt,
                   const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
  if(NucDeEx::Utils::fVerbose>0) std::cout << "NucDeExDeexcitationPholeDoDeex()  Z=" << Z << "   N=" << N << endl;
  EventID++;

  // --- Save event level info (except for shell & status) --- //
  EventInfo.InitParameters();
  SaveEventLevelInfo(Zt,Nt,Z,N,Ex,mom);

  // store target info
  Z_target   = Z;
  N_target   = N;
  mom_target = mom;
  nuc_target = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z_target,N_target);
  name_target = (string)nuc_target->name;

  double random = NucDeEx::Random::random();

  //-----------11B------------//
  if(Z==5 && N==6){
  //--------------------------//
    int index=0;
    double Br_integ=0;
    for(int i=0;i<Nlevel_p32_11B;i++){
      Br_integ += Br_p32_11B[i];
      if(random<Br_integ){
        index=i;
        break;
      }
    }
    if(index==0){ // g.s.
      AddGSNucleus(Z,N,mom);
    }else{ // excited state
      // set paremeters for boost calculation
      decay_mode=0; 
      Ex_target  = E_p32_11B[index];
      Ex_daughter=0;
      mass_particle=0;
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      mass_daughter=mass_target;
      name_daughter = name_target;
      Z_daughter = Z_target;
      N_daughter = N_target;
      nuc_daughter = nuc_target;
      S = nuc_target->S[decay_mode];
      Qvalue = Ex_target - S - Ex_daughter;
      Decay(1); // breakflag on
    }
  //-----------11C------------//
  }else if(Z==6 && N==5){
  //--------------------------//
    int index=0;
    double Br_integ=0;
    for(int i=0;i<Nlevel_p32_11C;i++){
      Br_integ += Br_p32_11C[i];
      if(random<Br_integ){
        index=i;
        break;
      }
    }
    if(index==0){ // g.s.
      AddGSNucleus(Z,N,mom);
    }else{ // excited state
      // set paremeters for boost calculation
      decay_mode=0; 
      Ex_target  = E_p32_11C[index];
      Ex_daughter=0;
      mass_particle=0;
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      mass_daughter=mass_target;
      name_daughter = name_target;
      Z_daughter = Z_target;
      N_daughter = N_target;
      nuc_daughter = nuc_target;
      S = nuc_target->S[decay_mode];
      Qvalue = Ex_target - S - Ex_daughter;
      Decay(1); // breakflag on
    }
  //-----------15N------------//
  }else if(Z==7 && N==8){
  //--------------------------//
    int index=0;
    double Br_integ=0;
    for(int i=0;i<Nlevel_p32_15N;i++){
      Br_integ += Br_p32_15N[i];
      if(random<Br_integ){
        index=i;
        break;
      }
    }
    if(index==0){ // gamma
      // set paremeters for boost calculation
      decay_mode=0; 
      Ex_target  = E_p32_15N[index];
      Ex_daughter=0;
      mass_particle=0;
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      mass_daughter=mass_target;
      name_daughter = name_target;
      Z_daughter = Z_target;
      N_daughter = N_target;
      nuc_daughter = nuc_target;
      S = nuc_target->S[decay_mode];
      Qvalue = Ex_target - S - Ex_daughter;
      Decay(1); // breakflag on
    }else if(index==1){ // gamma but multiple -> read talys
      EventInfo = deex_talys->DoDeex(Zt,Nt,Z,N,E_p32_15N[index],mom);
    }else if(index==2){
      // set paremeters for boost calculation
      decay_mode=2; // proton
      Ex_target  = E_p32_15N[index];
      Ex_daughter=0;
      mass_particle = fTDatabasePDG->GetParticle(NucDeEx::PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      Z_daughter = Z_target-1;
      N_daughter = N_target;
      mass_daughter = ElementMassInMeV(Z_daughter+N_daughter, Z_daughter);
      nuc_daughter = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z_daughter,N_daughter);
      name_daughter = nuc_daughter->name;
      S = nuc_target->S[decay_mode];
      Qvalue = Ex_target - S - Ex_daughter;
      Decay(1); // breakflag on
    }else{
      exit(1);
    }
  //-----------15O------------//
  }else if(Z==8 && N==7){
  //--------------------------//
    int index=0;
    double Br_integ=0;
    for(int i=0;i<Nlevel_p32_15O;i++){
      Br_integ += Br_p32_15O[i];
      if(random<Br_integ){
        index=i;
        break;
      }
    }
    if(index==0){ // gamma
      // set paremeters for boost calculation
      decay_mode=0; 
      Ex_target  = E_p32_15O[index];
      Ex_daughter=0;
      mass_particle=0;
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      mass_daughter=mass_target;
      name_daughter = name_target;
      Z_daughter = Z_target;
      N_daughter = N_target;
      S = nuc_target->S[decay_mode];
      Qvalue = Ex_target - S - Ex_daughter;
      Decay(1); // breakflag on
    }else if(index==1||index==2){
      // set paremeters for boost calculation
      decay_mode=2; // proton
      Ex_target  = E_p32_15O[index];
      Ex_daughter=0;
      mass_particle = fTDatabasePDG->GetParticle(NucDeEx::PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      Z_daughter = Z_target-1;
      N_daughter = N_target;
      mass_daughter = ElementMassInMeV(Z_daughter+N_daughter, Z_daughter);
      nuc_daughter = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z_daughter,N_daughter);
      name_daughter = nuc_daughter->name;
      S = nuc_target->S[decay_mode];
      Qvalue = Ex_target - S - Ex_daughter;
      Decay(1); // breakflag on
    }
  }else{ 
    exit(1);
  }
  return EventInfo;
}
