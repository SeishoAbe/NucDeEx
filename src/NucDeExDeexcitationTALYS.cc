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
#include "NucDeExDeexcitationTALYS.hh"

///////////////////////////
NucDeExDeexcitationTALYS::NucDeExDeexcitationTALYS(): ldmodel(2), parity_optmodall(1){};
///////////////////////////

///////////////////////////
NucDeExDeexcitationTALYS::NucDeExDeexcitationTALYS(const int ld, const bool p_o)
///////////////////////////
{
  ldmodel=ld;
  parity_optmodall=p_o;
  NucDeEx::Utils::SetPATH();
  NucDeEx::Utils::NucleusTable->ReadTables(0); // should be after SetPATH
}

#ifdef INCL_DEEXCITATION_NUCDEEX
///////////////////////////
NucDeExDeexcitationTALYS::NucDeExDeexcitationTALYS(const int ld, const bool p_o, G4INCL::Config *config)
///////////////////////////
{ 
  ldmodel=ld;
  parity_optmodall=p_o;
  NucDeEx::Utils::NucleusTable->ReadTables(0); // should be after SetPATH
}
#endif

/////////////////////////////////////////////
NucDeExEventInfo NucDeExDeexcitationTALYS::DoDeex(const int Zt, const int Nt,
                   const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
  if(NucDeEx::Utils::fVerbose>0) std::cout << "NucDeExDeexcitationTALYS::DoDeex()" << std::endl;
  EventID++;

  RESET:

  // --- Save event level info (except for shell & status) --- //
  EventInfo.InitParameters();
  SaveEventLevelInfo(Zt,Nt,Z,N,Ex,mom);
  EventInfo.fShell=1;

  // --- store target info --- //
  Z_target   = Z;
  N_target   = N;
  Ex_target  = Ex;
  if(Ex_target<0) Ex_target=0;
  mom_target = mom;
  nuc_target = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z_target,N_target);
  if(nuc_target==NULL){ // nothing can be done for this case...(no mass profile)
    if(NucDeEx::Utils::fVerbose>0){
      std::cout << "We don't have deexcitation profile for this nucleus: "
           << name_target.c_str() << std::endl;
    }
    EventInfo.fStatus=0;
    return EventInfo;
  }
  name_target = (string) nuc_target->name;

  // --- Read ROOT file --- //
  if( ! OpenROOT(Zt,Nt,Z,N) ){ // No ROOT file
    if(NucDeEx::Utils::fVerbose>0){
      std::cout << "We don't have deexcitation profile for this nucleus: " << name_target.c_str() << std::endl;
    }
    AddGSNucleus(Z,N,mom);
    EventInfo.fStatus=0;
    return EventInfo;
  }

  int status=1;
  
  // Loop until zero excitation energy or null nuc_daughter ptr
  // Use private members (parameters) named as "_target"
  while(true){// <- infinite loop. There is break point
    if(NucDeEx::Utils::fVerbose>0){
      std::cout << "### " << name_target << ",   Ex = " << Ex_target << "   mom_target: "; mom_target.Print();
    }

    // --- Get (TGraph*) br based on name_target
    if(GetBrTGraph(name_target)){ // TGraph found
      // --- Determine decay mode 
      //     Return: The same as array in NucDeExConsts.hh
      decay_mode = DecayMode(Ex_target); 

      // --- Get nearest Ex bin (TGraph point) and then get (TGraph*) br_ex
      if(GetBrExTGraph(name_target, Ex_target, decay_mode)==0){ // g.s.
        AddGSNucleus(Z,N,mom);
        return EventInfo;
      }

      // --- Determine daughter excitation energy
      int Ex_daughter_point;
      DaughterExPoint(&Ex_daughter,&Ex_daughter_point);

      // --- Get mass using ROOT libraries
      if(decay_mode<=2){ // obtained from TDatabasePDG
        mass_particle = NucDeEx::Utils::fTDatabasePDG->GetParticle(NucDeEx::PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
      }else{ // obtained from TGeoElementRN
        int a_particle = (NucDeEx::PDG_particle[decay_mode]%1000)/10;
        int z_particle = ((NucDeEx::PDG_particle[decay_mode]%1000000)-a_particle*10)/10000;
        mass_particle = ElementMassInMeV(a_particle,z_particle);
        if(NucDeEx::Utils::fVerbose>1) std::cout << "a_particle = " << a_particle << "   z_particle = " << z_particle << std::endl;
      }
      mass_target = ElementMassInMeV(Z_target+N_target, Z_target);
      mass_daughter = ElementMassInMeV(Z_daughter+N_daughter, Z_daughter);
      if(mass_target<0 || mass_daughter<0){// this rarely happens...
        if(NucDeEx::Utils::fVerbose>0){
          std::cout << "Cannot find " << name_daughter << " in TGeoElementRN" << std::endl;
          std::cout << "Call DoDeex() again!" << std::endl;
        }
        goto RESET; // call this fuction again
      }
    }else{ // no tgraph found -> gamma emission to g.s.
      if(NucDeEx::Utils::fVerbose>0){
        std::cout << "Cannot find TGraph" << std::endl;
        std::cout << "Force gamma decay" << std::endl;
      }
      decay_mode=0; 
      Ex_daughter=0;
      mass_particle=0;
      mass_daughter=mass_target;
    }

    // --- Get separation E and Qvalue 
    S = nuc_target->S[decay_mode];
    Qvalue = Ex_target - S - Ex_daughter;
    if(Qvalue<0) Qvalue=0;
    if(NucDeEx::Utils::fVerbose>1){
      std::cout << "S = " << S << std::endl;
      std::cout << "Qvalue = " << Qvalue << std::endl;
    }

    bool breakflag=0;
    if(nuc_daughter==NULL || Ex_daughter==0) breakflag=1;
    
    // --- Calculate kinematics --- //
    if(Decay(breakflag)==0) status=-1; 
      // if breakflag==0, it does not save daughter nucleus
      // return 0: Something suspicious thing happens in energy conservation

    if(breakflag) break;

    // end of while loop: daughter -> target 
    Z_target  = Z_daughter;
    N_target  = N_daughter;
    Ex_target = Ex_daughter;
    mass_target = mass_daughter;
    mom_target  += mom_daughter;
    nuc_target = nuc_daughter;
    name_target = (string)nuc_target->name;
    if(NucDeEx::Utils::fVerbose>0) std::cout << std::endl;

    // --- Need this to release memory of TGraph
    //    TGraph memory looks not relased only by closing & deleting root file...
    DeleteTGraphs();
  }



  rootf->Close();
  delete rootf;
  EventInfo.fStatus = status;
  return EventInfo;
}


/////////////////////////////////////////////
int NucDeExDeexcitationTALYS::DecayMode(const double Ex)
/////////////////////////////////////////////
{
  // --- Determine decay mode 
  // ---- Return: (int)decay_mode_r 

  double Br[NucDeEx::num_particle]={0};
  double Br_sum=0;
  for(int p=0;p<NucDeEx::num_particle;p++){
    Br[p] = g_br[p]->Eval(Ex);
    if(Br[p]<0) Br[p]=0;
    Br_sum += Br[p];
    if(NucDeEx::Utils::fVerbose>1){
      std::cout << "Br(" << NucDeEx::particle_name[p].substr(0,1) << ") = " << Br[p] << std::endl;
    }
  }

  double Br_integ=0;
  double random = NucDeEx::Random::random();
  int decay_mode_r=-1;
  for(int p=0;p<NucDeEx::num_particle;p++){
    Br[p] /= Br_sum; // Normalize Br just in case.  
    Br_integ += Br[p];
    if(NucDeEx::Utils::fVerbose>1){
      std::cout << "Br_integ(" << NucDeEx::particle_name[p].substr(0,1) << ") = " << Br_integ << std::endl;
    }
    if(Br_integ>random){
      decay_mode_r=p;
      break; 
    } 
  }
  if(decay_mode_r<0){
    std::cerr << "Unexpected decay_mode = " << decay_mode_r << std::endl;
    exit(1);
  }

  // --- Set Z and N for daughter 
  Z_daughter = Z_target;
  N_daughter = N_target;
  if(decay_mode_r==1){
    N_daughter--;
  }else if(decay_mode_r==2){
    Z_daughter--;
  }else if(decay_mode_r==3){
    Z_daughter--;
    N_daughter--;
  }else if(decay_mode_r==4){
    Z_daughter--;
    N_daughter-=2;
  }else if(decay_mode_r==5){
    Z_daughter-=2;
    N_daughter--;
  }else if(decay_mode_r==6){
    Z_daughter-=2;
    N_daughter-=2;
  }
  nuc_daughter = NucDeEx::Utils::NucleusTable->GetNucleusPtr(Z_daughter,N_daughter);
  if(nuc_daughter!=NULL) name_daughter = nuc_daughter->name;
  else{
    os.str("");
    os << Z_daughter+N_daughter << NucDeEx::Utils::NucleusTable->nuc_name[Z_daughter];
    name_daughter = os.str();
  }

  if(NucDeEx::Utils::fVerbose>1){ 
    std::cout << "DecayMode(): Random = " << random << " : " << name_target.c_str() << " --> " << NucDeEx::particle_name[decay_mode_r] << " + ";
    std::cout << name_daughter.c_str() << std::endl;
  }

  return decay_mode_r;
}

/////////////////////////////////////////////
bool NucDeExDeexcitationTALYS::OpenROOT(const int Zt,const int Nt, const int Z, const int N)
/////////////////////////////////////////////
{
  os.str("");
  os << NucDeEx::Utils::NUCDEEX_ROOT << "/output/";
  // single nucleon hole
  if( (Zt==Z && Nt==N+1) || (Zt==Z+1 && Nt==N) ){ 
    if(Zt==6&&Nt==6) os << "12C/";
    else if(Zt==8&&Nt==8) os << "16O/";
    else return 0; // not supported
  }
  os << "Br_" << name_target.c_str() << "_ldmodel" << ldmodel;
  if(parity_optmodall) os << "_parity_optmodall";
  os << ".root"; 
  //
  rootf = new TFile(os.str().c_str(),"READ");
  if(! rootf->IsOpen()) return 0;
  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "OpenRoot: " << os.str().c_str() << std::endl;
  }
  return 1;
}

/////////////////////////////////////////////
bool NucDeExDeexcitationTALYS::GetBrTGraph(const string st)
/////////////////////////////////////////////
{
  for(int p=0;p<NucDeEx::num_particle;p++){
    os.str("");
    os << "g_" << st.c_str() << "_br_" << p;
    g_br[p] = (TGraph*) rootf->Get(os.str().c_str());
    if(g_br[p]==0) return 0; // no tgraph
  }
  return 1;
}


/////////////////////////////////////////////
int NucDeExDeexcitationTALYS::GetBrExTGraph(const string st, const double ex_t, const int mode)
/////////////////////////////////////////////
{
  double ex,br;
  double diff_ex=0;
  int point=0;
  for(point=0;point<g_br[mode]->GetN();point++){
    g_br[mode]->GetPoint(point,ex,br);
    if(ex>ex_t) break;
    diff_ex = abs(ex-ex_t);
  }
  if(abs(ex-ex_t)>diff_ex) point--;
  if(point==g_br[mode]->GetN()) point--;
  if(point<0) point=0;
  os.str("");
  os << "g_" << st.c_str() << "_br_ex_" << mode << "_" << point;
  g_br_ex = (TGraph*) rootf->Get(os.str().c_str());
  if(g_br_ex==0) return -1; // no tgraph

  if(g_br_ex->GetN()==0){ //no point in the tgraph -> get next point
    point++;
    os.str("");
    os << "g_" << st.c_str() << "_br_ex_" << mode << "_" << point;
    g_br_ex = (TGraph*) rootf->Get(os.str().c_str());
    if(g_br_ex==0) return -1; // no tgraph
  }
  g_br[mode]->GetPoint(point,ex,br);

  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "GetBrExTGraph(): nearest_point = " << point << ",  ex at the point = " << ex
      << ", diff_Ex = " << abs(ex-ex_t) << ", diff_ex(previous) = " << diff_ex << std::endl;

  }
  if(ex==0 && point==0) return 0; // g.s.
  return 1;
}


/////////////////////////////////////////////
bool NucDeExDeexcitationTALYS::DaughterExPoint(double *d_Ex, int *d_point)
/////////////////////////////////////////////
{
  double ex=0, br=0;
  double Br_sum=0;
  for(int p=0;p<g_br_ex->GetN();p++){
    g_br_ex->GetPoint(p,ex,br);
    Br_sum += br;
  }
  double Br_integ=0;
  double random = NucDeEx::Random::random();
  int point=0;
  for(point=0;point<g_br_ex->GetN();point++){
    g_br_ex->GetPoint(point,ex,br);
    br /= Br_sum; // Normalize Br just in case.
    Br_integ += br;
    if(Br_integ>random) break;
  }
  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "DaughterExPoint: Random = " << random << " : ex = " << ex
         << ",   point = " << point << std::endl;
  }

  *d_Ex=ex;
  *d_point=point;

  return 1;
}


/////////////////////////////////////////////
void NucDeExDeexcitationTALYS::DeleteTGraphs()
/////////////////////////////////////////////
{
  for(int p=0;p<NucDeEx::num_particle;p++){
    delete g_br[p];
    g_br[p]=0;
  }
  delete g_br_ex;
  g_br_ex=0;
}
