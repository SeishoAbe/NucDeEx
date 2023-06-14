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
	eventID=0;
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
void Deexcitation::DoDeex(const int Z, const int N, const double Ex, const TVector3* mom)
/////////////////////////////////////////////
{
	cout << endl << "###################################" << endl;
	cout << "Deexcitation::DoDeex(" << Z << "," << N << "," 
			 << Ex << ")  eventID=" << eventID << endl;
	cout << "###################################" << endl;

	// --- Initialization --- //
	InitParticleVector();
	if(mom==0) mom_target.SetXYZ(0,0,0); // decay at rest
	else       mom_target.SetXYZ(mom->X(),mom->Y(),mom->Z());
	// store target info. we don't want to change original value.
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
	
	// Loop until zero excitation energy or null nuc_daughter ptr
	// Use private members (parameters) named as "_target"
	while(true){// <- infinite loop. There is break point
		cout << "### " << name_target << ",   Ex = " << Ex_target;
		cout << "     mom_target: "; mom_target.Print();

		// --- Get (TGraph*) br based on name_target
		GetBrTGraph(name_target); 
		
		// --- Determine decay mode 
		//     Return: The same as array in consts.hh
		//     It also sets Z_daughter, Z_daughter
		decay_mode = DecayMode(Ex_target); 

		// --- Get nearest Ex bin (TGraph point) and then get (TGraph*) br_ex
		GetBrExTGraph(name_target, Ex_target, decay_mode);

		// --- Determine daughter excitation energy
		int Ex_daughter_point;
		DaughterExPoint(&Ex_daughter,&Ex_daughter_point);

		// --- Get mass using ROOT libraries
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

		// --- Get separation E and Qvalue 
		S = nuc_target->S[decay_mode];
		Qvalue = Ex_target - S - Ex_daughter;
		if(Qvalue<0) Qvalue=0;
		if(verbose>0){
			cout << "S = " << S << endl;
			cout << "Qvalue = " << Qvalue << endl;
		}

		bool breakflag=0;
		if(nuc_daughter==NULL || Ex_daughter==0) breakflag=1;
		
		// --- Calculate kinematics --- //
		Decay(breakflag); // if breakflag==0, it does not save daughter nucleus

		if(breakflag) break;

		// end of while loop: daughter -> target 
		Z_target  = Z_daughter;
		N_target  = N_daughter;
		Ex_target = Ex_daughter;
		mass_target = mass_daughter;
		mom_target  += mom_daughter;
		//kE_target  = kE_daughter;
		nuc_target = _nucleus_table->GetNucleusPtr(Z_target,N_target);
		name_target = (string)nuc_target->name;
		cout << endl;
	}
		
	rootf->Close();
	delete rootf;
	eventID++;
}


/////////////////////////////////////////////
int Deexcitation::DecayMode(const double Ex)
/////////////////////////////////////////////
{
	// --- Determine decay mode 
	// ---- Return: (int)decay_mode 

	double Br[num_particle]={0};
	double Br_sum=0;
	for(int p=0;p<num_particle;p++){
		Br[p] = g_br[p]->Eval(Ex);
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
		cout << "DecayMode(): Random = " << random << " : " << name_target.c_str() << " --> " << particle_name[decay_mode] << " + ";
		if(nuc_daughter!=NULL) cout << nuc_daughter->name << endl;
		else cout << Z_daughter+N_daughter << _nucleus_table->nuc_name[Z_daughter] << endl;
	}

	return decay_mode;
}


/////////////////////////////////////////////
bool Deexcitation::DaughterExPoint(double *d_Ex, int *d_point)
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
		cout << "DaughterExPoint: Random = " << random << " : ex = " << ex
		     << ",   point = " << point << endl;
	}

	*d_Ex=ex;
	*d_point=point;

	return 1;
}

/////////////////////////////////////////////
void Deexcitation::Decay(bool breakflag)
/////////////////////////////////////////////
{
	cout << "Deexcitation::Decay()" << endl;

	cout << "mass_target   = " << mass_target << endl;
	cout << "mass_particle = " << mass_particle << endl;
	cout << "mass_daughter = " << mass_daughter << endl;
	
	// ---- CM frame --- //
	// --- Calculate kinematics (CM)
	//		mass used in the following calculation should include excitation energy
	double mass_ex_daughter = mass_daughter + Ex_daughter;

	double cmMomentum = std::sqrt(Qvalue*(Qvalue + 2.*mass_particle)*
															(Qvalue + 2.*mass_ex_daughter)*
															(Qvalue + 2.*mass_particle + 2.*mass_ex_daughter) )/
															(Qvalue + mass_particle + mass_ex_daughter)/2.;
	double kE_particle = sqrt( pow(cmMomentum,2) + pow(mass_particle,2) ) - mass_particle;
	double kE_daughter = sqrt( pow(cmMomentum,2) + pow(mass_ex_daughter,2) ) - mass_ex_daughter;
	double kE_sum = kE_particle + kE_daughter; // just for check
	if( Qvalue>0 && (kE_sum - Qvalue)/Qvalue > check_criteria){
		cerr << "Error @ DeexcitationD::Decay: Something wroing in kinematics calculation" << endl;
		abort();
	}
	if(verbose>0){
		cout << "cmMomentum  = " << cmMomentum << endl;
		cout << "kE_particle = " << kE_particle << endl;
		cout << "kE_daughter = " << kE_daughter << endl;
		cout << "kE_sum      = " << kE_sum << " <- consistent with Qvalue = " << Qvalue << endl;
	}
	
	// --- Detemine momentum direction (uniform) (CM)
	// --- Then set momentum vectors
	double costheta = 2.*rndm->Rndm()-1.; // [-1, 1]
	double sintheta = sqrt( 1. - pow(costheta,2) );
	double phi      = 2*TMath::Pi()*rndm->Rndm(); // [0,2pi]
	TVector3 dir( sintheta*cos(phi), sintheta*sin(phi), costheta );
	mom_particle = -1*cmMomentum*dir; // -P
	mom_daughter = cmMomentum*dir; // P
	if(verbose>0){
		cout << "dir:          "; dir.Print();
		cout << "mom_particle: ";mom_particle.Print();
		cout << "mom_daughter: ";mom_daughter.Print();
	}
	

	// ---- CM -> LAB --- //
	// --- Get total energy of target (i.e., CM frame) in the LAB frame for boost
	//     This is "Total CM energy" (in the LAB frame)
	//		 This should not include excitation energy
	double totalE_target = sqrt( pow(mom_target.Mag(),2) + pow(mass_target,2) );

	// --- Store info as Particle. then boost it
	Particle p_particle(PDG_particle[decay_mode],
											mass_particle,
											mom_particle,
											verbose);
	double totalE_particle_bef = p_particle.totalE();
	p_particle.Boost(totalE_target,mom_target);// BOOST!
	double totalE_particle_aft = p_particle.totalE();

	Particle p_daughter(PDGion(Z_daughter,N_daughter),
											mass_daughter, // w/o excitation E
											mom_daughter,
											verbose);
	double totalE_daughter_bef = p_daughter.totalE();
	p_daughter.Boost(totalE_target,mom_target);
	double totalE_daughter_aft = p_daughter.totalE();

	// 
	double kE_target = totalE_target - mass_target;
	double totalE_ex_target = totalE_target + Ex_target; // w/ excitation E
	cout << "totalE_target = " << totalE_target << endl;
	cout << "kE_target = " << kE_target << endl;
	cout << "Ex_target = " << Ex_target << endl;
	cout << "totalE_ex_target = " << totalE_ex_target << endl;

	double totalE_bef = totalE_particle_bef + totalE_daughter_bef;
	double totalE_aft = totalE_particle_aft + totalE_daughter_aft;
	cout << "totalE_bef = " << totalE_bef << endl;
	cout << "totalE_aft = " << totalE_aft << endl;
	cout << "  diff = " << totalE_aft-totalE_bef << endl;

	double totalE_ex_bef = totalE_bef + Ex_daughter; // w/ excitation E
	double totalE_ex_aft = totalE_aft + Ex_daughter;
	cout << "totalE_ex_bef = " << totalE_ex_bef << endl;
	cout << "totalE_ex_aft = " << totalE_ex_aft << endl;

	// --- Check Energy conservation (4dim energy including momentum)
	//       Fundamental energy conservation 
	//		   (Total energy in LAB) = (Total energy in CM after boost)
	// i.e., (Total energy of target w/ ex in LAB) 
	//         = (Total energy of particle after boost) + (Total energy of daughter w/ ex after boost)
	//						The last two terms are calculated from CM at first, and then boosted.
	if(totalE_ex_target>0 && (totalE_ex_aft-totalE_ex_target)/totalE_ex_target>check_criteria){
		cerr << "ERROR: @ Deexcitation:Decay(): Energy is not conserved..." << endl;
		abort();
	}
	_particle.push_back(p_particle);

	// DoDeex loop will be end -> Save daughter
	if(breakflag)_particle.push_back(p_daughter);
}

/////////////////////////////////////////////
double Deexcitation::ElementMassInMeV(TGeoElementRN* ele)
/////////////////////////////////////////////
{
	double mass = ele->MassNo()*(TGeoUnit::amu_c2/TGeoUnit::MeV)
									 + ele->MassEx(); // (MeV)
	double mass_amu = mass/(TGeoUnit::amu_c2/TGeoUnit::MeV);
	if(verbose>0){
		cout << "mass (MeV) = " << mass 
				 << "    (amu) = " << mass_amu << endl;
	}
	return mass;
}

/////////////////////////////////////////////
int Deexcitation::PDGion(int Z, int N)
/////////////////////////////////////////////
{
	int pdg= 1e9 + Z*1e4 + (Z+N)*1e1;
	return pdg;
}

/////////////////////////////////////////////
bool Deexcitation::OpenROOT(const char* name)
/////////////////////////////////////////////
{
	os.str("");
	os << getenv("TALYS_WORK") << "/output/Br_" << name << ".root"; 
	rootf = new TFile(os.str().c_str(),"READ");
	if(! rootf->IsOpen()) return 0;
	return 1;
}

/////////////////////////////////////////////
bool Deexcitation::GetBrTGraph(const string st)
/////////////////////////////////////////////
{
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "g_" << st.c_str() << "_br_" << p;
		g_br[p] = (TGraph*) rootf->Get(os.str().c_str());
	}
	return 1;
}

/////////////////////////////////////////////
int Deexcitation::GetBrExTGraph(const string st, const double ex_t, const int mode)
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
	if(verbose>0){
		cout << "GetBrExTGraph(): nearest_point = " << point << ",  diff_Ex = " << abs(ex-ex_t) << ", diff_ex(previous) = " << diff_ex << endl;
	}
	os.str("");
	os << "g_" << st.c_str() << "_br_ex_" << mode << "_" << point;
	g_br_ex = (TGraph*) rootf->Get(os.str().c_str());

	return point;
}

/////////////////////////////////////////////
void Deexcitation::InitParticleVector()
/////////////////////////////////////////////
{
	_particle.clear();
	vector<Particle>().swap(_particle); // memory release
}
