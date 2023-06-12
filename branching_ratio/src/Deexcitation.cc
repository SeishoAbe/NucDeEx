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
#include <TParticlePDG.h>
#include <TParticle.h>
#include <TGeoPhysicalConstants.h>

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
	rndm = new TRandom3(0); // 0: seed is from time
	pdg = new TDatabasePDG();
	geo = new TGeoManager("test","test");
	element_table = geo->GetElementTable();
}

///////////////////////////
Deexcitation::~Deexcitation()
///////////////////////////
{
	delete _nucleus_table;
	delete rndm;
	delete pdg;
	delete geo;
}

/////////////////////////////////////////////
void Deexcitation::DoDeex(int Z,int N, double Ex)
/////////////////////////////////////////////
{
	cout << endl << "###################################" << endl;
	cout << "Deexcitation::DoDeex(" << Z << "," << N << "," 
			 << Ex << ")" << endl;
	cout << "###################################" << endl;
	
	// store target info. we don't want to change original value...
	Z_target   = Z;
	N_target   = N;
	Ex_target  = Ex;
	nuc_target = _nucleus_table->GetNucleusPtr(Z_target,N_target);
	name_target = (string)nuc_target->name;

	// Read ROOT file
	if( ! OpenROOT(name_target.c_str()) ){
		cout << "We don't have deexcitation profile for this nucleus: " 
				 << name_target.c_str() << endl;
		return;
	}
	
	// loop until zero excitation energy
	// Use parameters named as "_target"
	while(Ex_target>0){
		cout << endl << "### " << name_target << ",   Ex = " << Ex_target << endl;
		GetBrTGraph(); // get TGraph based on name_target
		
		// --- Determine decay mode 
		//     Return: The same as array in consts.hh
		//     It also sets Z_daughter, Z_daughter
		int decay_mode = DecayMode(Ex_target); 

		// --- Get nearest Ex bin (TGraph point) & get (TGraph*) br_ex
		int Ex_target_point = NearestExPoint(Ex_target,decay_mode); 

		// --- Determine daughter excitation energy
		int Ex_daughter_point;
		DaughterExPoint(Ex_daughter,Ex_daughter_point);

		// get mass info using ROOT libraries
		if(decay_mode<=2){ // obtained from TDatabasePDG
			mass_particle = pdg->GetParticle(PDG_particle[decay_mode])->Mass()/TGeoUnit::MeV;//GeV2MeV
		}else{ // obtained from TGeoElementRN
			int a_particle = (PDG_particle[decay_mode]%1000)/10;
			int z_particle = ((PDG_particle[decay_mode]%1000000)-a_particle*10)/10000;
			TGeoElementRN* element_rn =	element_table->GetElementRN(a_particle,z_particle); // (A,Z)
			cout << "a_particle = " << a_particle << "   z_particle = " << z_particle << endl;
			mass_particle = ElementMassInMeV(element_rn);
		}
		mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
		mass_daughter = ElementMassInMeV(element_table->GetElementRN(Z_daughter+N_daughter, Z_daughter));

		// get separation E & Qvalue 
		S = nuc_target->S[decay_mode];
		Qvalue = Ex_target - S - Ex_daughter;
		if(Qvalue<0) Qvalue=0;
		if(verbose>0){
			cout << "S = " << S << endl;
			cout << "Qvalue = " << Qvalue << endl;
		}

		Decay();

		// end of while loop: daughter -> target 
		if(nuc_daughter==NULL) break;
		Ex_target = Ex_daughter;
		Z_target  = Z_daughter;
		N_target  = N_daughter;
		nuc_target = _nucleus_table->GetNucleusPtr(Z_target,N_target);
		name_target = (string)nuc_target->name;
	}
		


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

	// --- Set Z and N for daughter 
	Z_daughter = Z_target;
	N_daughter = N_target;
	if(decay_mode==1){
		N_daughter--;
	}else if(decay_mode==2){
		Z_daughter--;
	}else if(decay_mode==3){
		Z_daughter--;
		N_daughter--;
	}else if(decay_mode==4){
		Z_daughter--;
		N_daughter-=2;
	}else if(decay_mode==5){
		Z_daughter-=2;
		N_daughter--;
	}else if(decay_mode==6){
		Z_daughter-=2;
		N_daughter-=2;
	}
	nuc_daughter = _nucleus_table->GetNucleusPtr(Z_daughter,N_daughter);

	if(verbose>0){ 
		cout << "Random = " << random << " : " << name_target.c_str() << " --> " << particle_name[decay_mode] << " + ";
		if(nuc_daughter!=NULL) cout << nuc_daughter->name << endl;
		else cout << Z_daughter+N_daughter << nuc_name[Z_daughter] << endl;
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
		cout << "Random = " << random << " : ex = " << ex
		     << ",   point = " << point << endl;
	}

	d_Ex=ex;
	d_point=point;

	return 1;
}

/////////////////////////////////////////////
vector<TParticle*> Deexcitation::Decay()
/////////////////////////////////////////////
{
	cout << "Deexcitation::Decay()" << endl;

	cout << "mass_target (MeV)   = " << mass_target << endl;
	cout << "mass_particle (MeV) = " << mass_particle << endl;
	cout << "mass_daughter (MeV) = " << mass_daughter << endl;
	// HOGEHOGE //


	vector<TParticle*> particle;
	return particle;
}

/////////////////////////////////////////////
double Deexcitation::ElementMassInMeV(TGeoElementRN* ele)
/////////////////////////////////////////////
{
	double mass = ele->MassNo()*(TGeoUnit::amu_c2/TGeoUnit::MeV)
									 + ele->MassEx();
	double mass_amu = mass/(TGeoUnit::amu_c2/TGeoUnit::MeV);
	if(verbose>0){
		cout << "mass (MeV) = " << mass 
				 << "    (amu) = " << mass_amu << endl;
	}
	return mass;
}


/////////////////////////////////////////////
const char* Deexcitation::PDGion(int Z, int N)
/////////////////////////////////////////////
{
	os.str("");
	os << "10" << setw(3) << setfill('0') << Z
						 << setw(3) << setfill('0') << Z+N << "0";
	return os.str().c_str();
}



/////////////////////////////////////////////
bool Deexcitation::OpenROOT(const char* name)
/////////////////////////////////////////////
{
	os.str("");
	os << getenv("TALYS_WORK_BR") << "/output/Br_" << name << ".root"; 
	rootf = new TFile(os.str().c_str(),"READ");
	if(! rootf->IsOpen()) return 0;
	return 1;
}

/////////////////////////////////////////////
bool Deexcitation::GetBrTGraph()
/////////////////////////////////////////////
{
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "g_" << name_target.c_str() << "_br_" << p;
		g_br[p] = (TGraph*) rootf->Get(os.str().c_str());
	}
	return 1;
}
