#include <iostream>
#include <fstream>

#include "vcworkC.h"
#include "posinnucC.h"
#include "vcvrtxC.h"
#include "necardC.h"
#include "neworkC.h"

#include <TVector3.h>
#include <TFile.h>
#include <TTree.h>

#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "NucDeExEventInfo.hh"
#include "NucDeExDeexcitation.hh"

NucDeExDeexcitation* nucdeex;
static bool nucdeex_initialized=false;

TFile* tmp;
TTree* tree;
float Ex;
int mode;
int Z, N;

void NucDeExInitialize(void);
float ExcitationEnergy(int Zt, int Nt);
TVector3 NucleusMomentum();
TVector3 mom, mom2;
//ofstream ofs_tmp;

extern "C"{
  int nucdeex_();
  // 0: Not supported
  // 1: OK (free p is included here)
  // -1: Error in NucDeEx
}

///////////////////////////
int nucdeex_()
///////////////////////////
{
  // --- 
  // free proton  -> Nothing to do
  if(posinnuc_.ibound==0) return 1;
  // currently QE only
  if(! ( abs(nework_.modene)==1 ||
       (abs(nework_.modene)>=51&&  abs(nework_.modene)<=52) ) ) return 0;
  // target should be nucleon
  if(vcwork_.ipvc[1]!=2112 && vcwork_.ipvc[1]!=2212) return 0;
  int Zt = neuttarget_.numbndp;
  int Nt = neuttarget_.numbndn;
  // 12C and 16O target nuclei only
  if(! ( (Zt==6 && Nt==6) || (Zt==8 && Nt==8) ) ) return 0;

  // --- determine Z and N of residual nucleus --- //
  Z=Zt, N=Nt;
  bool flag_CC=0;// 1: CC, 0:NC
  if(abs(nework_.modene)<30) flag_CC=1;
  bool flag_nu=0;// 1: neutrino, 0: antineutrino
  if(vcwork_.ipvc[0]>0) flag_nu=1;
  //
  if(flag_CC){
    if(flag_nu){  // nu + n -> l- + p
      Z++;
      N--;
    }else{ // nubar + p -> l+ + n
      Z--;
      N++;
    }
  }
  // FSI 
  for(int i=3;i<vcwork_.nvc;i++){
    if(NucDeEx::Utils::fVerbose>=1){
      cout << "### " << i << endl;
      cout << vcwork_.ipvc[i] << "   icrnvc=" << vcwork_.icrnvc[i] << endl;
      cout << "iflgvc=" << vcwork_.iflgvc[i] << "   iorgvc=" << vcwork_.iorgvc[i] << endl;
    }
    if(vcwork_.icrnvc[i]!=1) continue; // remove pre-FSI
    if(vcwork_.ipvc[i]==2112) N--;
    else if(vcwork_.ipvc[i]==2212) Z--;
  }
  Z = std::max(0,Z);
  N = std::max(0,N);
  if(NucDeEx::Utils::fVerbose>=1){
    cout << "Z = " << Z << "   N = " << N << endl;
  }

  // Calculate input variables to NucDeEx
  NucDeExInitialize();
  Ex = ExcitationEnergy(Zt,Nt);
  mode = nework_.modene;
  //tree->Fill();

  TVector3 mom_nucleus = NucleusMomentum(); // mom of residual nucleus
  //ofs_tmp << nework_.modene << " " << Ex << "  " << mom_nucleus.Mag() << endl;
  NucDeExEventInfo theNucDeExResult = nucdeex->DoDeex(Zt,Nt,Z,N,Ex,mom_nucleus);
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
    vcwork_.ivtfvc[vcwork_.nvc+n_deex]  = 1;  // final vertex
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
  }
  vcwork_.nvc += n_deex;

  if(NucDeEx::Utils::fVerbose>=2){
    std::cout << "Zt=" << Zt << " Nt=" << Nt << "  mom=" << mom_nucleus.Mag() << " Ex=" << Ex << endl;
    std::cout << "vcwork nvc " << vcwork_.nvc << std::endl;
    std::cout << "vcwork posvc(" << vcwork_.posvc[0] << "," << vcwork_.posvc[1] << "," << vcwork_.posvc[2] << ")" << std::endl;
    std::cout << "###" << std::endl;
    for(int i=0;i<vcwork_.nvc;i++){
      std::cout << "ipvc=" << vcwork_.ipvc[i] << " icrnvc=" << vcwork_.icrnvc[i] << " iflgvc=" << vcwork_.iflgvc[i] << endl;
      std::cout << "mass=" << vcwork_.amasvc[i] 
                << "  " << sqrt( pow(vcwork_.pvc[i][0],2) + pow(vcwork_.pvc[i][1],2) + pow(vcwork_.pvc[i][2],2) ) << endl;
      //std::cout << "vcwork posivc("
      //          << vcwork_.posivc[i][0] << ","
      //          << vcwork_.posivc[i][1] << ","
      //          << vcwork_.posivc[i][2] << ")" << endl;
      //std::cout << "vcwork posfvc("
      //          << vcwork_.posfvc[i][0] << ","
      //          << vcwork_.posfvc[i][1] << ","
      //          << vcwork_.posfvc[i][2] << ")" << endl;
      //std::cout << "posinnuc posnuc("
      //          << posinnuc_.posnuc[i][0] << ","
      //          << posinnuc_.posnuc[i][1] << ","
      //          << posinnuc_.posnuc[i][2] << ")" << endl;
    }
  }
  //for(int i=0;i<vcvrtx_.nvtxvc;i++){
  //  std::cout << "vcvrtx pvtxvc("
  //            << vcvrtx_.pvtxvc[i][0] << ","
  //            << vcvrtx_.pvtxvc[i][1] << ","
  //            << vcvrtx_.pvtxvc[i][2] << ")" << endl;
  //}
  
  return theNucDeExResult.fStatus;
}


///////////////////////////
float ExcitationEnergy(int Zt, int Nt)
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

  float S; // separation E
  if(Zt==6 && Nt==6) S = NucDeEx::Utils::NucleusTable->GetNucleusPtr("12C")->S[2];
  else if(Zt==8 && Nt==8) S = NucDeEx::Utils::NucleusTable->GetNucleusPtr("16O")->S[2];
  else abort();

  Ex = MissE-S;

  return Ex;
}

///////////////////////////
TVector3 NucleusMomentum()
///////////////////////////
{
  mom.SetXYZ(vcwork_.pvc[1][0],vcwork_.pvc[1][1],vcwork_.pvc[1][2]);
  if(abs(nework_.modene)==2){ // 2p2h
    mom2.SetXYZ(vcwork_.pvc[2][0],vcwork_.pvc[2][1],vcwork_.pvc[2][2]);
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

  //ofs_tmp.open("tmp.txt");

  nucdeex = new NucDeExDeexcitation();
  // random seed (temporary)
  NucDeEx::Random::SetSeed(1);
  // verosity
  int neut_quiet = neutcard_.quiet;
  if(neut_quiet==0) NucDeEx::Utils::fVerbose=2;
  else if(neut_quiet==1) NucDeEx::Utils::fVerbose=1;
  else if(neut_quiet==2) NucDeEx::Utils::fVerbose=0;

  if(NucDeEx::Utils::fVerbose>0) std::cout << "NucDeEx initialized" << std::endl;
  nucdeex_initialized=true;
/*
  tmp = new TFile("tmp.root","RECREATE");
  tree = new TTree("tree","tree");
  tree->Branch("Ex",&Ex,"Ex/F");
  tree->Branch("Z",&Z,"Z/I");
  tree->Branch("N",&N,"N/I");
  tree->Branch("mode",&mode,"mode/I");
  tree->SetAutoSave(1);
*/

  return;
}
