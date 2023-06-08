#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <algorithm>

#include "Deexcitation.hh"
#include "Nucleus.hh"
#include "consts.hh"
#include <TGraph.h>
#include <TRandom3.h>

using namespace std;

///////////////////////////
Deexcitation::Deexcitation()
///////////////////////////
{
	// Prepare nuc table
	_nucleus_table = new NucleusTable();
	if(!_nucleus_table->ReadTables()){
		cerr << "Fatal Error" << endl;
		abort();
	}
	//if(!ReadROOT()) abort();
	verbose=0;
	ran = new TRandom3(0);
}


///////////////////////////
void Deexcitation::DoDeex(int Z,int N, double Ex)
///////////////////////////
{
	Z_target=Z;
	N_target=N;
	name_target = (string)_nucleus_table->GetNucleusPtr(Z_target,N_target)->name;

	os.str("");
	os << getenv("TALYS_WORK_BR") << "/output/Br_" << name_target.c_str() << ".root";
	rootf = new TFile(os.str().c_str(),"READ");
	tree = (TTree*) rootf->Get("tree");
	float Br[num_particle]={0};
	float Br_sum=0;
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "g_" << name_target.c_str() << "_br_" << p;
		g_br[p] = (TGraph*) rootf->Get(os.str().c_str());
		Br[p] = g_br[p]->Eval(Ex);
		Br_sum += Br[p];
		if(verbose>0){
			cout << "Br(" << particle_name[p].substr(0,1) << ") = " << Br[p] << endl;
		}
	}

	// --- Determine 1st decay mode
	//  (Normalize Br just in case)
	float Br_rand=0;
	float random = ran->Rndm();
	int particle_decay;
	for(int p=0;p<num_particle;p++){
		Br[p] /= Br_sum;
		Br_rand += Br[p];
		if(verbose>0){
			cout << "Br_rand(" << particle_name[p].substr(0,1) << ") = " << Br_rand << endl;
		}
		if(Br_rand>random){
			particle_decay=p;
			break;
		}
	}
	if(verbose>0){
		cout << "Random = " << random << " --> " << particle_name[particle_decay] << endl;
	}
	

	// --- Get br_ex TGraph
	double d_ex=0;
	int point=0;
	double ex, br;
	for(point=0;point<g_br[particle_decay]->GetN();point++){
		g_br[particle_decay]->GetPoint(point,ex,br);
		if(ex>Ex) break;
		d_ex = abs(ex-Ex);
	}
	if(abs(ex-Ex)>d_ex) point--;
	if(verbose>1){
		cout << "point = " << point << "  " << abs(ex-Ex) << "  prev_de_ex = " << d_ex << endl;
	}
	os.str("");
	os << "g_" << name_target.c_str() << "_br_ex_" << particle_decay << "_" << point;
	g_br_ex = (TGraph*) rootf->Get(os.str().c_str());


	// Determine daughter level ex
	Br_sum=0;
	for(int p=0;p<g_br_ex->GetN();p++){
		g_br_ex->GetPoint(p,ex,br);
		Br_sum += br;
	}
	Br_rand=0;
	random = ran->Rndm();
	int daughter_point=0;
	for(int p=0;p<g_br_ex->GetN();p++){
		g_br_ex->GetPoint(p,ex,br);
		br /= Br_sum;
		Br_rand += br;
		if(Br_rand>random){
			daughter_point = p;
			break;
		}
	}
	if(verbose>0){
		cout << "Random = " << random << " Goes to " << daughter_point << "  ex = " << ex << endl;
	}





	rootf->Close();
	delete rootf;
}




///////////////////////////
bool Deexcitation::ReadROOT()
///////////////////////////
{
	os.str("");
	os << getenv("TALYS_WORK_BR") << "/output/Br_11B.root";  // FIXME tentative
	rootf = new TFile(os.str().c_str(),"READ");
	tree = (TTree*) rootf->Get("tree");
	/*
	tree->SetBranchAddress("name",&name);
	tree->SetBranchAddress("Z",&Z);
	tree->SetBranchAddress("N",&N);
	tree->SetBranchAddress("Ex_bin",&Ex_bin);
	for(int i=0;i<tree->GetEntries();i++){ // nuc loop
		tree->GetEntry(i);
		cout << name->c_str() << endl;
		*/
		/*
		for(int p=0;p<num_particle;p++){
			os.str("");
			os << "g_" << name.c_str() << "_br_" << p;
			g_br[i][p] = (TGraph*) rootf->Get(os.str().c_str());
			for(int bin=0;bin<Ex_bin;bin++){
				os.str("");
				os << "g_" << name.c_str() << "_br_ex_" << p << "_"<< bin;
				g_br_ex[i][p][bin] = (TGraph*) rootf->Get(os.str().c_str());
			}
		}
		*/
	//}

	return 1;
}
