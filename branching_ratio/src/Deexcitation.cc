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
	verbose=0;
	rndm = new TRandom3(0);
}

///////////////////////////
Deexcitation::~Deexcitation()
///////////////////////////
{
	delete rndm;
	delete _nucleus_table;
}

/////////////////////////////////////////////
void Deexcitation::DoDeex(int Z,int N, double Ex)
/////////////////////////////////////////////
{
	cout << "Deexcitation::DoDeex(" << Z << "," << N << "," << Ex << ")" << endl;
	
	// store target info. we don't want to change original value...
	Z_target=Z;
	N_target=N;
	Ex_target=Ex;
	name_target = (string)_nucleus_table->GetNucleusPtr(Z_target,N_target)->name;

	// Read ROOT file
	if( ! ReadROOT(name_target.c_str()) ){
		cout << "We don't have deexcitation profile for this nucleus: " 
				 << name_target.c_str() << endl;
		return;
	}

	//while(Ex_target>0){
		
		// --- Determine decay mode 
		int decay_mode = DecayMode(Ex_target); 

		// --- Get nearest Ex bin (TGraph point) & get (TGraph*) br_ex
		int nearest_ex_point = NearestExPoint(Ex_target,decay_mode); 

		// --- Determine daughter excitation energy
		int daughter_ex_point;
		double daughter_ex;
		DaughterExPoint(daughter_ex,daughter_ex_point);


	//}

	rootf->Close();
	delete rootf;
}


/////////////////////////////////////////////
int Deexcitation::DecayMode(double Ex)
/////////////////////////////////////////////
{
	// --- Determine decay mode 
	// ---- Return: (int)decay_mode 

	double Br[num_particle]={0};
	double Br_sum=0;
	for(int p=0;p<num_particle;p++){
		Br[p] = g_br[p]->Eval(Ex_target);
		Br_sum += Br[p];
		if(verbose>0){
			cout << "Br(" << particle_name[p].substr(0,1) << ") = " << Br[p] << endl;
		}
	}

	double Br_integ=0;
	double random = rndm->Rndm();
	int decay_mode;
	for(int p=0;p<num_particle;p++){
		Br[p] /= Br_sum; // Normalize Br just in case.
		Br_integ += Br[p];
		if(verbose>0){
			cout << "Br_integ(" << particle_name[p].substr(0,1) << ") = " << Br_integ << endl;
		}
		if(Br_integ>random){
			decay_mode=p;
			break; 
		} 
	}
	if(verbose>0){
		cout << "Random = " << random << " --> " << particle_name[decay_mode] << endl;
	}

	return decay_mode;
}


/////////////////////////////////////////////
int Deexcitation::NearestExPoint(double Ex, int decay_mode)
/////////////////////////////////////////////
{
	double ex,br;
	double diff_ex=0;
	int point=0;
	for(point=0;point<g_br[decay_mode]->GetN();point++){
		g_br[decay_mode]->GetPoint(point,ex,br);
		if(ex>Ex_target) break;
		diff_ex = abs(ex-Ex_target);
	}
	if(abs(ex-Ex_target)>diff_ex) point--;
	if(verbose>0){
		cout << "nearest_point = " << point << ",  diff_Ex = " << abs(ex-Ex_target) << ", diff_ex(previous) = " << diff_ex << endl;
	}
	os.str("");
	os << "g_" << name_target.c_str() << "_br_ex_" << decay_mode << "_" << point;
	g_br_ex = (TGraph*) rootf->Get(os.str().c_str());

	return point;
}


/////////////////////////////////////////////
bool Deexcitation::DaughterExPoint(double &d_Ex, int &d_point)
/////////////////////////////////////////////
{
	double ex, br;
	double Br_sum=0;
	for(int p=0;p<g_br_ex->GetN();p++){
		g_br_ex->GetPoint(p,ex,br);
		Br_sum += br;
	}
	double Br_integ=0;
	double random = rndm->Rndm();
	int point=0;
	for(point=0;point<g_br_ex->GetN();point++){
		g_br_ex->GetPoint(point,ex,br);
		br /= Br_sum; // Normalize Br just in case.
		Br_integ += br;
		if(Br_integ>random) break;
	}
	if(verbose>0){
		cout << "Random = " << random << " --> ex = " << ex
		     << ",   point = " << point << endl;
	}

	d_Ex=ex;
	d_point=point;

	return 1;
}














/////////////////////////////////////////////
bool Deexcitation::ReadROOT(const char* name)
/////////////////////////////////////////////
{
	os.str("");
	os << getenv("TALYS_WORK_BR") << "/output/Br_" << name << ".root"; 
	rootf = new TFile(os.str().c_str(),"READ");
	if(! rootf->IsOpen()) return 0;
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "g_" << name << "_br_" << p;
		g_br[p] = (TGraph*) rootf->Get(os.str().c_str());
	}
	/*
	tree = (TTree*) rootf->Get("tree");
	tree->SetBranchAddress("name",&name);
	tree->SetBranchAddress("Z",&Z);
	tree->SetBranchAddress("N",&N);
	tree->SetBranchAddress("Ex_bin",&Ex_bin);
	*/
	return 1;
}
