#include <iostream>

#include "vcworkC.h"
#include "posinnucC.h"
#include "vcvrtxC.h"
#include "necardC.h"
#include "neworkC.h"

#include <TVector3.h>

#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "NucDeExEventInfo.hh"
#include "NucDeExDeexcitation.hh"

NucDeExDeexcitation* nucdeex;
static bool nucdeex_initialized=false;

void NucDeExInitialize(void);
float ExcitationEnergy();
TVector3 NucleusMomentum();

extern "C"{
  int nucdeex_();
}

///////////////////////////
int nucdeex_()
///////////////////////////
{
  //std::cout << vcwork_.ipvc[1] << "   ibound=" << posinnuc_.ibound << std::endl;
  if(posinnuc_.ibound==0) return 0; 
  if(vcwork_.ipvc[1]!=2112 && vcwork_.ipvc[1]!=2212) return 0;

  NucDeExInitialize();

  float Ex = ExcitationEnergy();
  TVector3 mom_nucleus = NucleusMomentum(); // mom of residual nucleus
  
  int Z = neuttarget_.numbndp;
  int N = neuttarget_.numbndn;

  NucDeExEventInfo theNucDeExResult = nucdeex->DoDeex(8,8,7,8,Ex,mom_nucleus);
  std::vector<NucDeExParticle> ParticleVector = theNucDeExResult.ParticleVector;
  int size=ParticleVector.size();
  int n_deex=0;
  for(int i=0;i<size;i++){
    NucDeExParticle particle = ParticleVector.at(i);
    if(!particle._flag) continue; // skip intermediate states
    // --- save parameter --- //
    vcwork_.ipvc[vcwork_.nvc+n_deex]    = particle._PDG;
    vcwork_.iorgvc[vcwork_.nvc+n_deex]  = 2;  // parent = initial nucleon
    vcwork_.iflgvc[vcwork_.nvc+n_deex]  = 10; // flag = deexcitation 
    vcwork_.icrnvc[vcwork_.nvc+n_deex]  = 1;  // flag = chase
    vcwork_.ivtivc[vcwork_.nvc+n_deex]  = 1;  // initial vertex
    vcwork_.ivtfvc[vcwork_.nvc+n_deex]  = 1; // final vertex
    vcwork_.timvc[vcwork_.nvc+n_deex]   = 0.; // starting time
    vcwork_.amasvc[vcwork_.nvc+n_deex]  = (float) particle._mass;
    vcwork_.pvc[vcwork_.nvc+n_deex][0]  = (float) particle._momentum.X();
    vcwork_.pvc[vcwork_.nvc+n_deex][1]  = (float) particle._momentum.Y();
    vcwork_.pvc[vcwork_.nvc+n_deex][2]  = (float) particle._momentum.Z();
    vcwork_.posivc[vcwork_.nvc+n_deex][0] = vcwork_.posvc[0];
    vcwork_.posivc[vcwork_.nvc+n_deex][1] = vcwork_.posvc[1];
    vcwork_.posivc[vcwork_.nvc+n_deex][2] = vcwork_.posvc[2];
    vcwork_.posfvc[vcwork_.nvc+n_deex][0] = vcwork_.posvc[0];
    vcwork_.posfvc[vcwork_.nvc+n_deex][1] = vcwork_.posvc[1];
    vcwork_.posfvc[vcwork_.nvc+n_deex][2] = vcwork_.posvc[2];
    //
    posinnuc_.posnuc[vcwork_.nvc+n_deex][0] = posinnuc_.posnuc[0][0];
    posinnuc_.posnuc[vcwork_.nvc+n_deex][1] = posinnuc_.posnuc[0][1];
    posinnuc_.posnuc[vcwork_.nvc+n_deex][2] = posinnuc_.posnuc[0][2];
    n_deex++;
    //std::cout << "vcwork ipvc " << vcwork_.ipvc[vcwork_.nvc+n_deex] << endl;
    //std::cout << "vcwork posivc("
    //          << vcwork_.posivc[vcwork_.nvc+n_deex][0] << ","
    //          << vcwork_.posivc[vcwork_.nvc+n_deex][1] << ","
    //          << vcwork_.posivc[vcwork_.nvc+n_deex][2] << ")" << endl;
    //std::cout << "vcwork posfvc("
    //          << vcwork_.posfvc[vcwork_.nvc+n_deex][0] << ","
    //          << vcwork_.posfvc[vcwork_.nvc+n_deex][1] << ","
    //          << vcwork_.posfvc[vcwork_.nvc+n_deex][2] << ")" << endl;
    //std::cout << "posinnuc posnuc("
    //          << posinnuc_.posnuc[vcwork_.nvc+n_deex][0] << ","
    //          << posinnuc_.posnuc[vcwork_.nvc+n_deex][1] << ","
    //          << posinnuc_.posnuc[vcwork_.nvc+n_deex][2] << ")" << endl;
  }
  vcwork_.nvc += n_deex;

  //std::cout << "################" << std::endl;
  //std::cout << "vcwork nvc " << vcwork_.nvc << std::endl;
  //std::cout << "vcwork posvc(" << vcwork_.posvc[0] << "," << vcwork_.posvc[1] << "," << vcwork_.posvc[2] << ")" << std::endl;
  //std::cout << "################" << std::endl;
  //for(int i=0;i<vcwork_.nvc;i++){
  //  std::cout << "vcwork ipvc " << vcwork_.ipvc[i] << "   iflgvc " << vcwork_.iflgvc[i] << endl;
  //  std::cout << "vcwork posivc("
  //            << vcwork_.posivc[i][0] << ","
  //            << vcwork_.posivc[i][1] << ","
  //            << vcwork_.posivc[i][2] << ")" << endl;
  //  std::cout << "vcwork posfvc("
  //            << vcwork_.posfvc[i][0] << ","
  //            << vcwork_.posfvc[i][1] << ","
  //            << vcwork_.posfvc[i][2] << ")" << endl;
  //  std::cout << "posinnuc posnuc("
  //            << posinnuc_.posnuc[i][0] << ","
  //            << posinnuc_.posnuc[i][1] << ","
  //            << posinnuc_.posnuc[i][2] << ")" << endl;
  //}
  //for(int i=0;i<vcvrtx_.nvtxvc;i++){
  //  std::cout << "vcvrtx pvtxvc("
  //            << vcvrtx_.pvtxvc[i][0] << ","
  //            << vcvrtx_.pvtxvc[i][1] << ","
  //            << vcvrtx_.pvtxvc[i][2] << ")" << endl;
  //}
  
  return theNucDeExResult.fStatus;
}


///////////////////////////
float ExcitationEnergy()
///////////////////////////
{
  float Ex,MissE;

  float Ev = sqrt( pow(vcwork_.pvc[0][0],2) + pow(vcwork_.pvc[0][1],2) + pow(vcwork_.pvc[0][2],2) );
  float mass_target;
  float mass_lepton, E_lepton;
  TVector3 mom_lepton(0,0,0);
  float mass_preFSI_nuc, E_preFSI_nuc;
  TVector3 mom_preFSI_nuc(0,0,0);
  if(abs(nework_.modene)!=2){  // non-2p2h
    mass_target = vcwork_.amasvc[1];
    mom_lepton.SetXYZ(vcwork_.pvc[2][0],vcwork_.pvc[2][1],vcwork_.pvc[2][2]);
    mass_lepton = vcwork_.amasvc[2];
    mass_preFSI_nuc = vcwork_.amasvc[3];
    mom_preFSI_nuc.SetXYZ(vcwork_.pvc[3][0],vcwork_.pvc[3][1],vcwork_.pvc[3][2]);
  }else{ // 2p2h
    mass_target = vcwork_.amasvc[1] + vcwork_.amasvc[2];
    mom_lepton.SetXYZ(vcwork_.pvc[3][0],vcwork_.pvc[3][1],vcwork_.pvc[3][2]);
    mass_lepton = vcwork_.amasvc[3];
    mass_preFSI_nuc = vcwork_.amasvc[4] + vcwork_.amasvc[5];
    mom_preFSI_nuc.SetXYZ(vcwork_.pvc[4][0],vcwork_.pvc[4][1],vcwork_.pvc[4][2]);
    TVector3 mom_preFSI_nuc2(vcwork_.pvc[5][0],vcwork_.pvc[5][1],vcwork_.pvc[5][2]);
    mom_preFSI_nuc += mom_preFSI_nuc2;
  }
  E_lepton = sqrt( pow(mass_lepton,2) + pow(mom_lepton.Mag(),2) );
  E_preFSI_nuc = sqrt( pow(mass_preFSI_nuc,2) + pow(mom_preFSI_nuc.Mag(),2) );

  MissE = Ev - E_lepton - E_preFSI_nuc + mass_target;

  Ex = MissE; // FIXME

  return Ex;
}

///////////////////////////
TVector3 NucleusMomentum()
///////////////////////////
{
  TVector3 mom(0,0,0);
  mom.SetXYZ(vcwork_.pvc[1][0],vcwork_.pvc[1][1],vcwork_.pvc[1][2]);
  if(abs(nework_.modene)==2){ // 2p2h
    TVector3 mom2(vcwork_.pvc[2][0],vcwork_.pvc[2][1],vcwork_.pvc[2][2]);
    mom += mom2;
  }
  // Assume the residual nucleus momentum is the exact opposite of one of target nucleon 
  // i.e., the target nucleus momentum is at rest.
  mom *= -1;

  return mom;
}

///////////////////////////
void NucDeExInitialize(void)
///////////////////////////
{
  if(nucdeex_initialized) return; // already initialized

  nucdeex = new NucDeExDeexcitation();
  NucDeEx::Random::SetSeed(1); // tmp
  
  // verosity
  int neut_quiet = neutcard_.quiet;
  if(neut_quiet==0) NucDeEx::Utils::fVerbose=2;
  else if(neut_quiet==1) NucDeEx::Utils::fVerbose=1;
  else if(neut_quiet==2) NucDeEx::Utils::fVerbose=0;
  nucdeex_initialized=true;

  if(NucDeEx::Utils::fVerbose>0) std::cout << "NucDeEx initialized" << std::endl;

  return;
}
