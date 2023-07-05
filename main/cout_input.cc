#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "consts.hh"
#include "NucleusTable.hh"

using namespace std;

int main(int argc, char* argv[]){
	if(argc<=2){
		cerr << argv[0] << " [Target nucleus]" << endl;
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

	ostringstream os;
	
  NucleusTable* nucleus_table = new NucleusTable();
  if(!nucleus_table->ReadTables()){
		cerr << "something wrong" << endl;
		return 1;
	}
	Nucleus* nuc_target = nucleus_table->GetNucleusPtr(argv[1]);
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
	ofstream ofs(os.str().c_str());

	ofs << "projectile 0" << endl;
	ofs << "element " << Z << endl;
	ofs << "mass " << A << endl;
	if( flag_jpi && (A==11 || A==15)) {
		ofs << "energy energy.1.2.p" << endl;
	}else{
		ofs << "energy energy" << endl;
	}
	ofs << endl;

	ofs << "bins 0" << endl;
	ofs << "ldmodel " << ldmodel << endl;
	if(parity) ofs << "parity y" << endl;
	else       ofs << "parity n" << endl;
	if(optmodall) ofs << "optmodall y" << endl;
	else          ofs << "optmodall n" << endl;
	ofs << endl;

	// maxlevelsbin
	ofs << "maxlevelstar " << nuc_target->maxlevelsbin << endl;
	for(int p=0;p<num_particle;p++){
		int Zt=Z, Nt=N;
		ofs << "maxlevelsbin " << particle_name[p].substr(0,1) << " ";
		if(particle_name[p].substr(0,1)=="n") Nt--;
		else if(particle_name[p].substr(0,1)=="p") Zt--;
		else if(particle_name[p].substr(0,1)=="d") Zt--,Nt--;
		else if(particle_name[p].substr(0,1)=="t") Zt--,Nt-=2;
		else if(particle_name[p].substr(0,1)=="h") Zt-=2,Nt--;
		else if(particle_name[p].substr(0,1)=="a") Zt-=2,Nt-=2;
		int At = Zt+Nt;
		int maxlevelsbin=0;
		for(int i=0;i<nucleus_table->GetNumofNuc();i++){
			Nucleus* nuc = nucleus_table->GetNucleusPtr(i);
			if(nuc->Z==Zt && nuc->N==Nt  && nuc->A==At){
				cout << nuc->name << endl;
				maxlevelsbin=nuc->maxlevelsbin;
			}
		}
		ofs << maxlevelsbin;
		ofs << endl;
	}
	ofs << "maxlevelsres 0" << endl;
	ofs << "ejectiles g n p a d t h" << endl;

	// output
	ofs << endl;
	ofs << "outdiscrete y" << endl;
	ofs << "outpopulation y" << endl;
	ofs << "outdecay y" << endl;
	ofs << "outlevels y" << endl;

	ofs.close();

	return 0;
}
