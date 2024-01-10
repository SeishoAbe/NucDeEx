#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "NucDeExConsts.hh"
#include "NucDeExNucleusTable.hh"

int main(int argc, char* argv[]){
  if(argc<=2){
    std::cerr << argv[0] << " [Target nucleus]" << std::endl;
    return 0;
  }
  int ldmodel=1; // default: 1 (constant temperature)
  bool parity=0, optmodall=0;
  bool flag_jpi=0; 
  // 0: jpi undefined energy
  // 1: jpi defined energy (only for single nucleon hole)
  if(argc>=3) ldmodel = atoi(argv[2]);
  if(argc>=4) parity  = (bool) atoi(argv[3]);
  if(argc>=5) optmodall = (bool) atoi(argv[4]);
  if(argc>=6) flag_jpi=(bool) atoi(argv[5]);

  std::ostringstream os;
  
  NucDeExNucleusTable* nucleus_table = new NucDeExNucleusTable();
  if(!nucleus_table->ReadTables()){
    std::cerr << "something wrong" << std::endl;
    return 1;
  }
  NucDeExNucleus* nuc_target = nucleus_table->GetNucleusPtr(argv[1]);
  const int Z=nuc_target->Z;
  const int N=nuc_target->N;
  const int A=Z+N;

  os.str("");
  os << "input/";
  if(flag_jpi){
    if(A==11) os << "12C";
    else if(A==15) os << "16O";
    else abort();
  }
  os << "/input_" << argv[1] << "_ldmodel" << ldmodel;
  if(parity) os << "_parity";
  if(optmodall) os << "_optmodall";
  std::ofstream ofs(os.str().c_str());

  ofs << "projectile 0" << std::endl;
  ofs << "element " << Z << std::endl;
  ofs << "mass " << A << std::endl;
  if( flag_jpi && (A==11 || A==15)) {
    ofs << "energy energy.1.2.p" << std::endl;
  }else{
    ofs << "energy energy" << std::endl;
  }
  ofs << std::endl;

  ofs << "bins 0" << std::endl;
  ofs << "ldmodel " << ldmodel << std::endl;
  if(parity) ofs << "parity y" << std::endl;
  else       ofs << "parity n" << std::endl;
  if(optmodall) ofs << "optmodall y" << std::endl;
  else          ofs << "optmodall n" << std::endl;
  ofs << std::endl;

  // maxlevelsbin
  ofs << "maxlevelstar " << nuc_target->maxlevelsbin << std::endl;
  for(int p=0;p<NucDeEx::num_particle;p++){
    int Zt=Z, Nt=N;
    ofs << "maxlevelsbin " << NucDeEx::particle_name[p].substr(0,1) << " ";
    if(NucDeEx::particle_name[p].substr(0,1)=="n") Nt--;
    else if(NucDeEx::particle_name[p].substr(0,1)=="p") Zt--;
    else if(NucDeEx::particle_name[p].substr(0,1)=="d") Zt--,Nt--;
    else if(NucDeEx::particle_name[p].substr(0,1)=="t") Zt--,Nt-=2;
    else if(NucDeEx::particle_name[p].substr(0,1)=="h") Zt-=2,Nt--;
    else if(NucDeEx::particle_name[p].substr(0,1)=="a") Zt-=2,Nt-=2;
    int At = Zt+Nt;
    int maxlevelsbin=0;
    for(int i=0;i<nucleus_table->GetNumofNuc();i++){
      NucDeExNucleus* nuc = nucleus_table->GetNucleusPtr(i);
      if(nuc->Z==Zt && nuc->N==Nt  && nuc->A==At){
        std::cout << nuc->name << std::endl;
        maxlevelsbin=nuc->maxlevelsbin;
      }
    }
    ofs << maxlevelsbin;
    ofs << std::endl;
  }
  ofs << "maxlevelsres 0" << std::endl;
  ofs << "ejectiles g n p a d t h" << std::endl;

  // output
  ofs << std::endl;
  ofs << "outdiscrete y" << std::endl;
  ofs << "outpopulation y" << std::endl;
  ofs << "outdecay y" << std::endl;
  ofs << "outlevels y" << std::endl;

  ofs.close();

  return 0;
}
