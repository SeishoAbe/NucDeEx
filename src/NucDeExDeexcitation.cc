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
#include "NucDeExDeexcitation.hh"

///////////////////////////
NucDeExDeexcitation::NucDeExDeexcitation():ldmodel(2),parity_optmodall(1)
///////////////////////////
{
  NucDeEx::Utils::SetPATH();
  Init();
}

///////////////////////////
NucDeExDeexcitation::NucDeExDeexcitation(const int ld, const bool p_o): ldmodel(ld), parity_optmodall(p_o)
///////////////////////////
{
  NucDeEx::Utils::SetPATH();
  Init();
}

#ifdef INCL_DEEXCITATION_NUCDEEX
///////////////////////////
NucDeExDeexcitation::NucDeExDeexcitation(const int ld, const bool p_o, G4INCL::Config *config): ldmodel(ld), parity_optmodall(p_o)
///////////////////////////
{ 
  NucDeEx::Utils::SetPATH(config);
  Init();
}
#endif

///////////////////////////
NucDeExDeexcitation::~NucDeExDeexcitation()
///////////////////////////
{
  delete deex_talys;
  delete deex_phole;
}

///////////////////////////
void NucDeExDeexcitation::Init()
///////////////////////////
{
  NucDeEx::Utils::NucleusTable->ReadTables(0); // should be after SetPATH
  deex_talys = new NucDeExDeexcitationTALYS(ldmodel,parity_optmodall);
  deex_phole = new NucDeExDeexcitationPhole();
  deex_phole->SetPtrTALYS(deex_talys);
}

/////////////////////////////////////////////
NucDeExEventInfo NucDeExDeexcitation::DoDeex(const int Zt, const int Nt,
                          const int Z, const int N, const int shell, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
  if(! ((Zt==6 && Nt==6 )||(Zt==8 && Nt==8)) ){
    std::cerr << "ERROR: This tool does not support the target nucleus" << std::endl;
    exit(1);
  }
  if(NucDeEx::Utils::fVerbose>0){
    std::cout << std::endl << "###################################" << std::endl;
    std::cout << "NucDeExDeexcitation::DoDeex(" << Zt << "," << Nt<< ","  << Z << "," << N 
         << "," << shell << "," << Ex << ")  EventID=" << EventID << std::endl;
    std::cout << "###################################" << std::endl;
  }
  

  // --- Save event level info --- //
  EventInfo.InitParameters();
  SaveEventLevelInfo(Zt,Nt,Z,N,Ex,mom);
  fShell=-1;

  // --- Call sub functions according to shell and nucleus conditions- --//
  if( (Zt==Z && Nt==N+1) || (Zt==Z+1 && Nt==N) ){
    // --- Single nucleon disapperance
    // determine shell level of hole
    if(shell>0) fShell=shell;
    else if(shell==0) fShell=ExtoShell(Zt,Nt,Ex);
    else exit(1);
    //
    if(fShell==3){ // p1/2-hole. nothing to do
      if(NucDeEx::Utils::fVerbose>0) std::cout << "(p1/2)-hole -> g.s." <<std::endl;
      AddGSNucleus(Z,N,mom);
    }else if(fShell==2){ // p3/2-hole 
      EventInfo = deex_phole->DoDeex(Zt,Nt,Z,N,Ex,mom); 
    }else if(fShell==1){// s1/2-hole read TALYS data
      EventInfo = deex_talys->DoDeex(Zt,Nt,Z,N,Ex,mom);
    }
  }else if( (Zt+Nt>Z+N) || (Zt+Nt == Z+N) ){
    // multi nucleon disapperanace or no change in Atomic number (Coherent scattering etc.)
    fShell=1;
    EventInfo = deex_talys->DoDeex(Zt,Nt,Z,N,Ex,mom);
  }else{ // Zt and Nt are larger than Z and N. This can be happen in the use in Geant4. Do nothing.
    if(NucDeEx::Utils::fVerbose>0){
      std::cerr << "Waring: Unexpected target & residual nuclei: "
                << "Zt = " << Zt << "   Nt = " << Nt << "  Z = " << Z << "   N = " << N << std::endl;
    }
    fShell=1;
    EventInfo.fStatus=0; // 0 -> not supported
  }

  if(fShell<0 || fShell>3){
    std::cerr << "ERROR: Unexpected shell level: shell = " << fShell << std::endl;
    exit(1);
  }

  // event level info
  EventInfo.fShell = fShell;

  EventID++;
  return EventInfo;
}


/////////////////////////////////////////////
int NucDeExDeexcitation::ExtoShell(const int Zt, const int Nt, const double Ex)
/////////////////////////////////////////////
{
  if(Zt==6&&Nt==6){ // 12C
    if(Ex>Ex_12C_s12) return 1; // s1/2-hole
    else return 2; //p3/2-hole
  }else if(Zt==8&&Nt==8){ // 16O
    if(Ex>Ex_16O_s12) return 1;
    else if(Ex>Ex_16O_p32) return 2;
    else return 3;
  }else abort();
  return -1;
}

