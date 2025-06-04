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
NucDeExDeexcitationPhole::NucDeExDeexcitationPhole(): version(2)
///////////////////////////
{
  SetParameters();
}
///////////////////////////
NucDeExDeexcitationPhole::NucDeExDeexcitationPhole(int v): version(v)
///////////////////////////
{
  SetParameters();
}

/////////////////////////////////////////////
void NucDeExDeexcitationPhole::SetParameters()
/////////////////////////////////////////////
{
  if(version==1){ // old method until v1.3
    // 11B* (Panin et al., Phys. Lett. B 753 204-210. Experimental data)
    Nlevel_p32_11B = 3; 
    E_p32_11B  = new double[Nlevel_p32_11B]{0.,2.125,5.020};
    Br_p32_11B = new double[Nlevel_p32_11B]{0.82,0.10,0.08};
    // 11C* (assume analogy of 11B*)
    //      the same Br, but energy is different
    Nlevel_p32_11C = 3;
    E_p32_11C  = new double[Nlevel_p32_11C]{0.,2.000,4.804};
    Br_p32_11C = new double[Nlevel_p32_11C]{0.82,0.10,0.08};
    // 15N* (Ejili, Phys. Rev. C 58, 3)
    //     Multiple gamma from 9.93 MeV will be described from TALYS data
    //     Deexcitation from 10.70 is taken from Ejiri (100% proton)
    Nlevel_p32_15N = 3;
    E_p32_15N  = new double[Nlevel_p32_15N]{6.32,9.93,10.70};
    Br_p32_15N = new double[Nlevel_p32_15N]{0.872,0.064,0.064};
    // 15O* (assume analogy of 15N*)
    Nlevel_p32_15O=3;
    E_p32_15O  = new double[Nlevel_p32_15O]{6.18,9.61,10.48};
    Br_p32_15O = new double[Nlevel_p32_15O]{0.872,0.064,0.064}; // guess
    // Br_p32_15O = new double[Nlevel_p32_15O]{1.,0,0}; // original Ejiri's value

  }else if(version==2 || version==3){ // new method from v2.1 (**only below separation E**)
    // 11B* (Panin et al., Phys. Lett. B 753 204-210. Experimental data)
    Nlevel_p32_11B = 3; 
    E_p32_11B  = new double[Nlevel_p32_11B]{0.,2.125,5.020};
    Br_p32_11B = new double[Nlevel_p32_11B]{0.82,0.10,0.08};
    // 11C* (assume analogy of 11B*)
    //      the same Br, but energy is different
    Nlevel_p32_11C = 3;
    E_p32_11C  = new double[Nlevel_p32_11C]{0.,2.000,4.804};
    Br_p32_11C = new double[Nlevel_p32_11C]{0.82,0.10,0.08};
    // 15N* (Yosoi PhD thesis & M. Leuschener et al., Phys. Rev. C 49, 955)
    //   - Neglecting 8.31 MeV 1/2+ and 9.05 MeV 1/2+
    //   - 5.27 MeV 5/2+ and 5.30 MeV 1/2+ are combined 
    Nlevel_p32_15N = 3;
    E_p32_15N  = new double[Nlevel_p32_15N]{5.27,6.32,9.93};
    Br_p32_15N = new double[Nlevel_p32_15N]{0.09,0.87,0.04};
    // 15O* (assume analogy of 15N*)
    //    - 5.18 MeV 1/2+ and 5.24 MeV 5/2+ are combined
    Nlevel_p32_15O=2;
    E_p32_15O  = new double[Nlevel_p32_15O]{5.18,6.18};
    Br_p32_15O = new double[Nlevel_p32_15O]{0.09,0.91};
  }else exit(1);
}

/////////////////////////////////////////////
NucDeExEventInfo NucDeExDeexcitationPhole::DoDeex(const int Zt, const int Nt,
                   const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
  if(NucDeEx::Utils::fVerbose>0) std::cout << "NucDeExDeexcitationPhole::DoDeex()  Z=" << Z << "   N=" << N << endl;
  EventID++;

  // --- Save event level info --- //
  EventInfo.InitParameters();
  SaveEventLevelInfo(Zt,Nt,Z,N,Ex,mom);

  // store target info
  Z_target   = Z;
  N_target   = N;
  mom_target = mom;
  nuc_target = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z_target,N_target);
  name_target = (string)nuc_target->name;

  if(version==1)      DoDeex_v1(Zt,Nt,Z,N,Ex,mom);
  else if(version==2) DoDeex_v2(Zt,Nt,Z,N,Ex,mom);
  else if(version==3){
    // v3 is only applicable for 11C and 11B. For others, use v2 instead.
    if((Z==5 && N==6)||(Z==6 && N==5)){
      DoDeex_v3(Zt,Nt,Z,N,Ex,mom);
    }else{
      DoDeex_v2(Zt,Nt,Z,N,Ex,mom);
    }
  }
  else exit(1);

  return EventInfo;
}


/////////////////////////////////////////////
void NucDeExDeexcitationPhole::DoDeex_v3(const int Zt, const int Nt,
                   const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
  // set paremeters for boost calculation
  decay_mode=0; 
  Ex_target  = Ex;
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
  for(int i=0;i<EventInfo.ParticleVector.size();i++){
      NucDeExParticle p = EventInfo.ParticleVector.at(i);
      cout << p._PDG << " " << p.kE() << endl;
  }
  if(EventInfo.ParticleVector.size()==1) EventInfo.fShell=3; // g.s.
  else EventInfo.fShell=2; // gamma discrete
}


/////////////////////////////////////////////
void NucDeExDeexcitationPhole::DoDeex_v2(const int Zt, const int Nt,
                   const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
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
      EventInfo.fShell=3;
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
      EventInfo.fShell=2;
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
      EventInfo.fShell=3;
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
      EventInfo.fShell=2;
    }
  //-----------15N------------//
  }else if(Z==7 && N==8){
  //--------------------------//
    if(Ex<E_p32_15N[0]){
      AddGSNucleus(Z,N,mom); // below 1st excited state
      EventInfo.fShell=3;
      return;
    }
    int index=0;
    double Br_integ=0;
    for(int i=0;i<Nlevel_p32_15N;i++){
      Br_integ += Br_p32_15N[i];
      if(random<Br_integ){
        index=i;
        break;
      }
    }
    if(index<=1){ // gamma
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
      EventInfo.fShell=2;
    }else if(index==2){ // multiple gamma
      EventInfo = deex_talys->DoDeex(Zt,Nt,Z,N,E_p32_15N[index],mom);
      EventInfo.fShell=2; //overwrite
    }else{
      exit(1);
    }
  //-----------15O------------//
  }else if(Z==8 && N==7){
  //--------------------------//
    if(Ex<E_p32_15O[0]){
      AddGSNucleus(Z,N,mom); // below 1st excited state
      EventInfo.fShell=3;
      return;
    }
    int index=0;
    double Br_integ=0;
    for(int i=0;i<Nlevel_p32_15O;i++){
      Br_integ += Br_p32_15O[i];
      if(random<Br_integ){
        index=i;
        break;
      }
    }
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
    EventInfo.fShell=2;
  }else{ 
    exit(1);
  }
}


/////////////////////////////////////////////
void NucDeExDeexcitationPhole::DoDeex_v1(const int Zt, const int Nt,
                   const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
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
      mass_particle = NucDeEx::Utils::fTDatabasePDG->GetParticle(NucDeEx::PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
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
      mass_particle = NucDeEx::Utils::fTDatabasePDG->GetParticle(NucDeEx::PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
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
  EventInfo.fShell=2; // overwrite
}
