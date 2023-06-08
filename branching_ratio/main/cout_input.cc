#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "consts.hh"
#include "NucleusTable.hh"

using namespace std;

int main(int argc, char* argv[]){
	if(argc!=2){
		cerr << argv[0] << " [Target nucleus]" << endl;
		return 0;
	}

	ostringstream os;

  NucleusTable* nucleus_table = new NucleusTable();
  if(!nucleus_table->ReadTables()){
		cerr << "something wrong" << endl;
		return 1;
	}

	os.str("");
	os << "input/input_" << argv[1];
	Nucleus* nuc_target = nucleus_table->GetNucleusPtr(argv[1]);

	ofstream ofs(os.str().c_str());
	ofs << "projectile 0" << endl;
	ofs << "element " << nuc_target->Z << endl;
	ofs << "mass " << nuc_target->A << endl;
	if(nuc_target->A==11 || nuc_target->A==15){
		ofs << "energy energy.1.2.p.txt" << endl;
	}else{
		ofs << "energy energy" << endl;
	}
	ofs << endl;

	ofs << "bins 0" << endl;
	// maxlevelsbin
	ofs << "maxlevelstar " << nuc_target->maxlevelsbin << endl;


	for(int p=0;p<num_particle;p++){
		int Z=nuc_target->Z;
		int N=nuc_target->N;
		ofs << "maxlevelsbin " << particle_name[p].substr(0,1) << " ";
		if(particle_name[p].substr(0,1)=="n") N--;
		else if(particle_name[p].substr(0,1)=="p") Z--;
		else if(particle_name[p].substr(0,1)=="d") Z--,N--;
		else if(particle_name[p].substr(0,1)=="t") Z--,N-=2;
		else if(particle_name[p].substr(0,1)=="h") Z-=2,N--;
		else if(particle_name[p].substr(0,1)=="a") Z-=2,N-=2;
		int A = Z+N;
		Nucleus* nuc;
		for(int i=0;i<nucleus_table->GetNumofNuc();i++){
			nuc = nucleus_table->GetNucleusPtr(i);
			if(nuc->Z==Z && nuc->N==N  && nuc->A==A) break;
		}
		cout << nuc->name << endl;
		ofs << nuc->maxlevelsbin;
		ofs << endl;
	}
	ofs << "maxlevelsres 0" << endl;
	ofs << "ejectiles g n p a d t h" << endl;

	// output
	ofs << endl;
	ofs << "outdiscrete y" << endl;
	ofs << "outpopulation y" << endl;
	//ofs << "outspectra y" << endl;
	ofs << "outdecay y" << endl;

	ofs.close();

	return 0;
}
