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

#include <TKey.h>

#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "NucDeExDeexcitationTALYS.hh"

///////////////////////////
NucDeExDeexcitationTALYS::NucDeExDeexcitationTALYS(): ldmodel(2), parity_optmodall(1)
///////////////////////////
{
  NucDeEx::Utils::SetPATH();
  NucDeEx::Utils::NucleusTable->ReadTables(0); // should be after SetPATH
  GetAllTGraph();
}

///////////////////////////
NucDeExDeexcitationTALYS::NucDeExDeexcitationTALYS(const int ld, const bool p_o): ldmodel(ld), parity_optmodall(p_o)
///////////////////////////
{
  NucDeEx::Utils::SetPATH();
  NucDeEx::Utils::NucleusTable->ReadTables(0); // should be after SetPATH
  GetAllTGraph();
}

#ifdef INCL_DEEXCITATION_NUCDEEX
///////////////////////////
NucDeExDeexcitationTALYS::NucDeExDeexcitationTALYS(const int ld, const bool p_o, G4INCL::Config *config): ldmodel(ld), parity_optmodall(p_o)
///////////////////////////
{ 
  NucDeEx::Utils::SetPATH(config);
  NucDeEx::Utils::NucleusTable->ReadTables(0); // should be after SetPATH
  GetAllTGraph();
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

  // --- Save event level info --- //
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
  if( (Zt==Z && Nt==N+1) || (Zt==Z+1 && Nt==N) ) fRootID = getRootID((name_target+"_1h").c_str());
  else                                           fRootID = getRootID(name_target.c_str());
  if(fRootID<0){ // No 
    if(NucDeEx::Utils::fVerbose>0){
      std::cout << "We don't have deexcitation profile for this nucleus: " << name_target.c_str() << std::endl;
    }
    AddGSNucleus(Z,N,mom);
    EventInfo.fStatus=0;
    return EventInfo;
  }
    
  // Loop until zero excitation energy or null nuc_daughter ptr
  // Use private members (parameters) named as "_target"
  while(true){// <- infinite loop. There is break point
    if(NucDeEx::Utils::fVerbose>0){
      std::cout << "### " << name_target << ",   Ex = " << Ex_target << "   mom_target: "; mom_target.Print();
    }

    // --- Get (TGraph*) br based on name_target
    fNucleusID = getNucleusID(fRootID,name_target.c_str());
    if(fNucleusID>=0){ // TGraph found
      // --- Determine decay mode 
      //     Return: The same as array in NucDeExConsts.hh
      decay_mode = DecayMode(Ex_target); 

      // --- Get nearest Ex bin (TGraph point) and then get (TGraph*) br_ex
      fPoint = GetBrExTGraph(name_target, Ex_target, decay_mode);
      if(fPoint==0){ // g.s.
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
    Decay(breakflag);

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
  }

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
    Br[p] = g_br_all[fRootID][fNucleusID][p]->Eval(Ex);
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
int NucDeExDeexcitationTALYS::GetBrExTGraph(const string st, const double ex_t, const int mode)
/////////////////////////////////////////////
{
  double ex,br;
  double diff_ex=0;
  int point=0;
  for(point=0;point<g_br_all[fRootID][fNucleusID][mode]->GetN();point++){
    g_br_all[fRootID][fNucleusID][mode]->GetPoint(point,ex,br);
    if(ex>ex_t) break;
    diff_ex = abs(ex-ex_t);
  }
  if(abs(ex-ex_t)>diff_ex) point--;
  if(point==g_br_all[fRootID][fNucleusID][mode]->GetN()) point--;
  if(point<0) point=0;
  //no point in the tgraph -> get next point
  if(g_br_ex_all[fRootID][fNucleusID][mode][point]->GetN()==0) point++;

  g_br_all[fRootID][fNucleusID][mode]->GetPoint(point,ex,br);

  if(NucDeEx::Utils::fVerbose>1){
    std::cout << "GetBrExTGraph(): nearest_point = " << point << ",  ex at the point = " << ex
      << ", diff_Ex = " << abs(ex-ex_t) << ", diff_ex(previous) = " << diff_ex << std::endl;

  }
  if(ex==0 && point==0) return 0; // g.s.
  return point;
}


/////////////////////////////////////////////
bool NucDeExDeexcitationTALYS::DaughterExPoint(double *d_Ex, int *d_point)
/////////////////////////////////////////////
{
  double ex=0, br=0;
  double Br_sum=0;
  for(int p=0;p<g_br_ex_all[fRootID][fNucleusID][decay_mode][fPoint]->GetN();p++){
    g_br_ex_all[fRootID][fNucleusID][decay_mode][fPoint]->GetPoint(p,ex,br);
    Br_sum += br;
  }
  double Br_integ=0;
  double random = NucDeEx::Random::random();
  int point=0;
  for(point=0;point<g_br_ex_all[fRootID][fNucleusID][decay_mode][fPoint]->GetN();point++){
    g_br_ex_all[fRootID][fNucleusID][decay_mode][fPoint]->GetPoint(point,ex,br);
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
void NucDeExDeexcitationTALYS::GetAllTGraph()
/////////////////////////////////////////////
{
  // --- Get list of Root files 
  char buf[256];
  int numfile=0;
  int pos;
  os.str("");
  os << "ls " << NucDeEx::Utils::NUCDEEX_ROOT << "/output/*.root";
  FILE* pipe = popen(os.str().c_str(),"r");
  std::string file[NucDeEx::bins];
  bool flag_1hole[NucDeEx::bins]={0};
  if(!pipe) exit(1);
  try{
    while(fgets(buf,sizeof buf,pipe) != NULL){
      file[numfile] = buf;
      pos= file[numfile].find("\n");
      if(pos!=std::string::npos){
        file[numfile] = file[numfile].substr(0,pos);
      }
      numfile++;
    }
  } catch (...){
    pclose(pipe);
    exit(1);
  }
  pclose(pipe);
  // 1hole
  os.str("");
  os << "ls " << NucDeEx::Utils::NUCDEEX_ROOT << "/output/*/*ldmodel2_parity*.root";
  pipe = popen(os.str().c_str(),"r");
  if(!pipe) exit(1);
  try{
    while(fgets(buf,sizeof buf,pipe) != NULL){
      file[numfile] = buf;
      pos= file[numfile].find("\n");
      if(pos!=std::string::npos){
        file[numfile] = file[numfile].substr(0,pos);
      }
      flag_1hole[numfile]=1; // flag on
      numfile++;
    }
  } catch (...){
    pclose(pipe);
    exit(1);
  }
  pclose(pipe);
  

  // --- Read Root files
  for(int i=0;i<numfile;i++){
    TFile* rootf = new TFile(file[i].c_str(),"READ");
    if(NucDeEx::Utils::fVerbose>0) std::cout << file[i] << std::endl;
    if(! rootf->IsOpen() ) exit(1);
    int pos1 = file[i].find("Br_");
    int pos2 = file[i].find("_ldmodel");
    if(pos1!=std::string::npos && pos2!=std::string::npos){
      std::string nucleus_name = file[i].substr(pos1+3,pos2-pos1-3);
      if(flag_1hole[i]) nucleus_name += "_1h";
      if(NucDeEx::Utils::fVerbose>2) std::cout << " -> " << nucleus_name << " (" << i << ")" << std::endl;
      // --- register nucleus of root file
      map_root.insert(std::make_pair(nucleus_name.c_str(),i));
    }else exit(1);
    TIter next(rootf->GetListOfKeys());
    TKey* key;
    int nucleus_id=0, nucleus_id_r=0; // _r is for counter.
    while( (key = (TKey*) next()) ){
      if(strstr(key->GetClassName(),"TGraph")){
        string name_obj = key->GetName();
        // --- remove pop tgraph
        pos = name_obj.find("pop");
        if(pos!=std::string::npos) continue;

        // --- now br or br_ex tgraph are selected. Set/Get index for [target]
        std::string name_target = name_obj.substr(2,name_obj.find("_br")-2);
        nucleus_id = getNucleusID(i,name_target.c_str());
        if(nucleus_id<0){
          map_nucleus[i].insert(std::make_pair(name_target,nucleus_id_r));
          nucleus_id = nucleus_id_r;
          nucleus_id_r++;
        }

        pos = name_obj.find("ex_");
        if(pos==std::string::npos){ // br tgraph
          // --- register if this target nucleus is not registered yet
          int index_particle = stoi(name_obj.substr(name_obj.length()-1));
          g_br_all[i][nucleus_id][index_particle] = (TGraph*) rootf->Get(name_obj.c_str());
          if(NucDeEx::Utils::fVerbose>2){
            std::cout << name_obj << " -> (" << i << ") " << name_target << " (" << nucleus_id << "), "
                      << index_particle << std::endl;
          }
        }else{ // br_ex tgraph
          int index_particle = stoi(name_obj.substr(pos+3,1));
          int index_ex_bin   = stoi(name_obj.substr(pos+5));
          g_br_ex_all[i][nucleus_id][index_particle][index_ex_bin] = (TGraph*) rootf->Get(name_obj.c_str());
          if(NucDeEx::Utils::fVerbose>2){
            std::cout << name_obj << " -> (" << i << ") "<< name_target << " (" << nucleus_id << "), "
                      << index_particle << ", " << index_ex_bin << std::endl;
          }
        }
      }
    } // end of key loop
    rootf->Close();
    delete rootf;
  } // end of root file loop
}


/////////////////////////////////////////////
int NucDeExDeexcitationTALYS::getRootID(const char* name)
/////////////////////////////////////////////
{
  it = map_root.find(name);
  if(it != map_root.end()){
    return (int) it->second;
  }else{
    return (int) -1;
  }
}

/////////////////////////////////////////////
int NucDeExDeexcitationTALYS::getNucleusID(int i, const char* name)
/////////////////////////////////////////////
{
  it = map_nucleus[i].find(name);
  if(it != map_nucleus[i].end()){
    return (int) it->second;
  }else{
    return -1;
  }
}
