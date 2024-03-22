#include <iostream>

#include "vcworkC.h"
#include "posinnucC.h"
#include "vcvrtxC.h"

#include "NucDeExUtils.hh"
#include "NucDeExDeexcitation.hh"

NucDeExDeexcitation* nucdeex;
static bool nucdeex_initialized=false;

void NucDeExInitialize(void);

extern "C"{
  int nucdeex_();
}


///////////////////////////
int nucdeex_()
t ///////////////////////////
{
  //std::cout << vcwork_.ipvc[1] << "   " << posinnuc_.ibound << std::endl;
  if(vcwork_.ipvc[1]!=2112 && vcwork_.ipvc[1]!=2212) return 0;
  if(posinnuc_.ibound==0) return 0; 

  NucDeExInitialize();
  TVector3 mom(0,0,0);
  int status = nucdeex->DoDeex(8,8,7,8,0,30,mom);
  int shell = nucdeex->GetShell();
  std::vector<NucDeExParticle> *particle = nucdeex->GetParticleVector();
  int size=particle->size();
  int n_deex=0;
  for(int i=0;i<size;i++){
    NucDeExParticle p = particle->at(i);
    if(!p._flag) continue;  // intermediate state 
    // --- save parameter --- //
    //p._momentum.Print();
    vcwork_.ipvc[vcwork_.nvc+n_deex]    = p._PDG;
    vcwork_.iorgvc[vcwork_.nvc+n_deex]  = 2;  // parent = initial nucleon
    vcwork_.iflgvc[vcwork_.nvc+n_deex]  = 10; // flag = deexcitation 
    vcwork_.icrnvc[vcwork_.nvc+n_deex]  = 1;  // flag = chase
    vcwork_.ivtivc[vcwork_.nvc+n_deex]  = 1;  // initial vertex
    vcwork_.ivtfvc[vcwork_.nvc+n_deex]  = 1; // final vertex
    vcwork_.timvc[vcwork_.nvc+n_deex]   = 0.; // starting time
    vcwork_.amasvc[vcwork_.nvc+n_deex]  = (float) p._mass;
    vcwork_.pvc[vcwork_.nvc+n_deex][0]  = (float) p._momentum.X();
    vcwork_.pvc[vcwork_.nvc+n_deex][1]  = (float) p._momentum.Y();
    vcwork_.pvc[vcwork_.nvc+n_deex][2]  = (float) p._momentum.Z();
    vcwork_.posivc[vcwork_.nvc+n_deex][0] = vcwork_.posvc[0];
    vcwork_.posivc[vcwork_.nvc+n_deex][1] = vcwork_.posvc[1];
    vcwork_.posivc[vcwork_.nvc+n_deex][2] = vcwork_.posvc[2];
    vcwork_.posfvc[vcwork_.nvc+n_deex][0] = vcwork_.posvc[0];
    vcwork_.posfvc[vcwork_.nvc+n_deex][1] = vcwork_.posvc[1];
    vcwork_.posfvc[vcwork_.nvc+n_deex][2] = vcwork_.posvc[2];
    // tmp
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
  return status;
}

///////////////////////////
void NucDeExInitialize(void)
///////////////////////////
{
  if(!nucdeex_initialized){
    std::cout << "NucDeEx initialized" << std::endl;
    nucdeex = new NucDeExDeexcitation(2,1);
    NucDeExUtils::SetSeed(1);
    NucDeExUtils::SetVerbose(1);
  }
  nucdeex_initialized=true;
  return;
}
