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

using namespace std;

///////////////////////////
Deexcitation::Deexcitation(const int ld, const bool p_o)
///////////////////////////
{
	// Prepare nuc table
	_nucleus_table = new NucleusTable();
	if(!_nucleus_table->ReadTables(0)){
		cerr << "Fatal Error" << endl;
		abort();
	}
	verbose=0;
	rndm = new TRandom3(0); // 0: seed is from time
	pdg = new TDatabasePDG();
	geo = new TGeoManager("test","test");
	element_table = geo->GetElementTable();
	eventID=0;
	ldmodel=ld;
	parity_optmodall=p_o;
	tree=0;
	for(int p=0;p<num_particle;p++){
		g_br[p]=0;
	}
	g_br_ex=0;
	_particle=0;
}

///////////////////////////
Deexcitation::~Deexcitation()
///////////////////////////
{
	delete _nucleus_table;
	delete rndm;
	delete pdg;
	delete geo;
	if(_particle!=0){
		_particle->clear();
		delete _particle;
	}
}

/////////////////////////////////////////////
int Deexcitation::DoDeex(const int Zt, const int Nt,
													const int Z, const int N, const int shell, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
	cout << endl << "###################################" << endl;
	cout << "Deexcitation::DoDeex(" << Zt << "," << Nt<< ","  << Z << "," << N 
			 << "," << shell << "," << Ex << ")  eventID=" << eventID << endl;
	cout << "###################################" << endl;
	eventID++;

	int status=-1;

	if(shell>0) _shell=shell;
	else if(shell==0) _shell=ExtoShell(Zt,Nt,Ex);
	else abort();

	if(! ((Zt==6 && Nt==6 )||(Zt==8 && Nt==8)) ){
		cerr << "This tool does not support the target nucleus" << endl;
		status=0;
	}

	// --- Call sub functions according to shell and nucleus conditions- --//
	if(Zt+Nt == Z+N || Ex<=0){
		// --- No change in nucleus (Coherent scattering etc.) or no nucleon emission
		//     Currently not supported. Nothing to do.
		// --- Negative Ex
		InitParticleVector();
		AddGSNucleus(Z,N,mom);
		status=0;
	}else if( (Zt==Z && Nt==N+1) || (Zt==Z+1 && Nt==N) ){
		// --- Single nucleon disapperance
		if(_shell==3){ 
			// p1/2-hole. nothing to do
			cout << "(p1/2)-hole" <<endl;
			InitParticleVector();
			AddGSNucleus(Z,N,mom);
			status=1;
		}else if(_shell==2){
			// p3/2-hole 
			status=DoDeex_p32(Zt,Nt,Z,N,mom); 
		}else if(_shell==1){
			// s1/2-hole read TALYS data
			status=DoDeex_talys(Zt,Nt,Z,N,Ex,mom);
		}else{
			cerr << "ERROR: Unexpected shell level: shell = " << _shell << endl;
			status=-1;
		}
	}else if(Zt+Nt>Z+N){
		// --- Multi-nucleon disapperance
		status=DoDeex_talys(Zt,Nt,Z,N,Ex,mom);
	}else{
		cerr << "ERROR: Unexpected target & residual nuclei" << endl;
		status=-1;
	}
	return status;
}

/////////////////////////////////////////////
int Deexcitation::DoDeex_talys(const int Zt, const int Nt,
													     const int Z, const int N, const double Ex, const TVector3& mom)
/////////////////////////////////////////////
{
	cout << "DoDeex_talys()" << endl;
	RESET:

	// --- Initialization --- //
	InitParticleVector();

	// store target info. we don't want to change original value.
	Z_target   = Z;
	N_target   = N;
	Ex_target  = Ex;
	mom_target = mom;
	nuc_target = _nucleus_table->GetNucleusPtr(Z_target,N_target);
	if(nuc_target==NULL){
		cout << "We don't have deexcitation profile for this nucleus: "
				 << name_target.c_str() << endl;
		//AddGSNucleus(Z,N,mom); // nothing can be done for this case...(no mass profile)
		return 0;
	}
	name_target = (string)nuc_target->name;

	// Read ROOT file
	if( ! OpenROOT(Zt,Nt,Z,N,0) ){
		cout << "We don't have deexcitation profile for this nucleus: " 
				 << name_target.c_str() << endl;
		AddGSNucleus(Z,N,mom);
		return 0;
	}
	
	// Loop until zero excitation energy or null nuc_daughter ptr
	// Use private members (parameters) named as "_target"
	while(true){// <- infinite loop. There is break point
		cout << "### " << name_target << ",   Ex = " << Ex_target;
		cout << "     mom_target: "; mom_target.Print();

		// --- Get (TGraph*) br based on name_target
		if(GetBrTGraph(name_target)){ // TGraph found
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
				mass_particle = pdg->GetParticle(PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
			}else{ // obtained from TGeoElementRN
				int a_particle = (PDG_particle[decay_mode]%1000)/10;
				int z_particle = ((PDG_particle[decay_mode]%1000000)-a_particle*10)/10000;
				TGeoElementRN* element_rn =	element_table->GetElementRN(a_particle,z_particle); // (A,Z)
				if(verbose>0) cout << "a_particle = " << a_particle << "   z_particle = " << z_particle << endl;
				mass_particle = ElementMassInMeV(element_rn);
			}
			// this rarely happens...
			if(element_table->GetElementRN(Z_daughter+N_daughter,Z_daughter) == NULL){
				cout << "Cannot find " << name_daughter << " in TGeoElementRN" << endl;
				cout << "Call DoDeex() again!" << endl;
				goto RESET; // call this fuc again
			}
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			mass_daughter = ElementMassInMeV(element_table->GetElementRN(Z_daughter+N_daughter, Z_daughter));
		}else{ // no tgraph found -> gamma emission to g.s.
			cout << "Cannot find TGraph" << endl;
			cout << "Force gamma decay" << endl;
			decay_mode=0; 
			Ex_daughter=0;
			mass_particle=0;
			mass_daughter=mass_target;
		}

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
		nuc_target = nuc_daughter;
		name_target = (string)nuc_target->name;
		cout << endl;

		// --- Need this to release memory of TGraph
		//		TGraph memory looks not relased only by closing & deleting root file...
		DeleteTGraphs();
	}

	rootf->Close();
	delete rootf;
	tree=0;

	return 1;
}

/////////////////////////////////////////////
int Deexcitation::DoDeex_p32(const int Zt, const int Nt,
													   const int Z, const int N, const TVector3& mom)
/////////////////////////////////////////////
{
	// --- Initialization --- //
	InitParticleVector();

	// store target info. we don't want to change original value.
	Z_target   = Z;
	N_target   = N;
	mom_target = mom;
	nuc_target = _nucleus_table->GetNucleusPtr(Z_target,N_target);
	name_target = (string)nuc_target->name;
	cout << "DoDeex_p32()" << endl;
	cout << "### " << name_target;
	cout << "     mom_target: "; mom_target.Print();

	double random = rndm->Rndm();

	//-----------11B------------//
	if(Z==5 && N==6){
	//--------------------------//
		int index=0;
		double Br_integ=0;
		for(int i=0;i<Nlevel_p32_11B;i++){
			Br_integ += Br_p32_11B[i];
			if(random<Br_integ){
				index=i;
				break;
			}
		}
		if(index==0){ // g.s.
			AddGSNucleus(Z,N,mom);
		}else{ // excited state
			// set paremeters for boost calculation
			decay_mode=0; 
			Ex_target  = E_p32_11B[index];
			Ex_daughter=0;
			mass_particle=0;
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			mass_daughter=mass_target;
			name_daughter = name_target;
			Z_daughter = Z_target;
			N_daughter = N_target;
			nuc_daughter = nuc_target;
			S = nuc_target->S[decay_mode];
			Qvalue = Ex_target - S - Ex_daughter;
			Decay(1); // breakflag on
		}
	//-----------11C------------//
	}else if(Z==6 && N==5){
	//--------------------------//
		int index=0;
		double Br_integ=0;
		for(int i=0;i<Nlevel_p32_11C;i++){
			Br_integ += Br_p32_11C[i];
			if(random<Br_integ){
				index=i;
				break;
			}
		}
		if(index==0){ // g.s.
			AddGSNucleus(Z,N,mom);
		}else{ // excited state
			// set paremeters for boost calculation
			decay_mode=0; 
			Ex_target  = E_p32_11C[index];
			Ex_daughter=0;
			mass_particle=0;
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			mass_daughter=mass_target;
			name_daughter = name_target;
			Z_daughter = Z_target;
			N_daughter = N_target;
			nuc_daughter = nuc_target;
			S = nuc_target->S[decay_mode];
			Qvalue = Ex_target - S - Ex_daughter;
			Decay(1); // breakflag on
		}
	//-----------15N------------//
	}else if(Z==7 && N==8){
	//--------------------------//
		int index=0;
		double Br_integ=0;
		for(int i=0;i<Nlevel_p32_15N;i++){
			Br_integ += Br_p32_15N[i];
			if(random<Br_integ){
				index=i;
				break;
			}
		}
		if(index==0){ // gamma
			// set paremeters for boost calculation
			decay_mode=0; 
			Ex_target  = E_p32_15N[index];
			Ex_daughter=0;
			mass_particle=0;
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			mass_daughter=mass_target;
			name_daughter = name_target;
			Z_daughter = Z_target;
			N_daughter = N_target;
			nuc_daughter = nuc_target;
			S = nuc_target->S[decay_mode];
			Qvalue = Ex_target - S - Ex_daughter;
			Decay(1); // breakflag on
		}else if(index==1){ // gamma but multiple -> read talys
			DoDeex_talys(Zt,Nt,Z,N,E_p32_15N[index],mom);
		}else if(index==2){
			// set paremeters for boost calculation
			decay_mode=2; // proton
			Ex_target  = E_p32_15N[index];
			Ex_daughter=0;
		  mass_particle = pdg->GetParticle(PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			Z_daughter = Z_target-1;
			N_daughter = N_target;
			mass_daughter = ElementMassInMeV(element_table->GetElementRN(Z_daughter+N_daughter, Z_daughter));
			nuc_daughter = _nucleus_table->GetNucleusPtr(Z_daughter,N_daughter);
			name_daughter = nuc_daughter->name;
			S = nuc_target->S[decay_mode];
			Qvalue = Ex_target - S - Ex_daughter;
			Decay(1); // breakflag on
		}else{
			abort();
		// return -1;
		}
	//-----------15O------------//
	}else if(Z==8 && N==7){
	//--------------------------//
		int index=0;
		double Br_integ=0;
		for(int i=0;i<Nlevel_p32_15O;i++){
			Br_integ += Br_p32_15O[i];
			if(random<Br_integ){
				index=i;
				break;
			}
		}
		if(index==0){ // gamma
			// set paremeters for boost calculation
			decay_mode=0; 
			Ex_target  = E_p32_15O[index];
			Ex_daughter=0;
			mass_particle=0;
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			mass_daughter=mass_target;
			name_daughter = name_target;
			Z_daughter = Z_target;
			N_daughter = N_target;
			S = nuc_target->S[decay_mode];
			Qvalue = Ex_target - S - Ex_daughter;
			Decay(1); // breakflag on
		}else if(index==1||index==2){
			// set paremeters for boost calculation
			decay_mode=2; // proton
			Ex_target  = E_p32_15O[index];
			Ex_daughter=0;
		  mass_particle = pdg->GetParticle(PDG_particle[decay_mode])->Mass()*1e3;// GeV2MeV
			mass_target = ElementMassInMeV(element_table->GetElementRN(Z_target+N_target, Z_target));
			Z_daughter = Z_target-1;
			N_daughter = N_target;
			mass_daughter = ElementMassInMeV(element_table->GetElementRN(Z_daughter+N_daughter, Z_daughter));
			nuc_daughter = _nucleus_table->GetNucleusPtr(Z_daughter,N_daughter);
			name_daughter = nuc_daughter->name;
			S = nuc_target->S[decay_mode];
			Qvalue = Ex_target - S - Ex_daughter;
			Decay(1); // breakflag on
		}
	}else{ 
		abort();
		// return -1;
	}
	return 1;
}

/////////////////////////////////////////////
void Deexcitation::AddGSNucleus(const int Z,const int N, const TVector3& mom)
/////////////////////////////////////////////
{
	mass_target = ElementMassInMeV(element_table->GetElementRN(Z+N, Z));
	nuc_target = _nucleus_table->GetNucleusPtr(Z,N);
	if(nuc_target==NULL || mass_target<0) return; // do nothing
	name_target = (string)nuc_target->name;
	Particle nucleus(PDGion(Z,N),
								   mass_target, // w/o excitation E
									 mom,
								   name_target, 
									 1,0,verbose); // track flag on // zero ex
	_particle->push_back(nucleus);
	cout << "AddGSNucleus(): " << name_target << endl;
}


/////////////////////////////////////////////
int Deexcitation::ExtoShell(const int Zt, const int Nt, const double Ex)
/////////////////////////////////////////////
{
	if(Zt==6&&Nt==6){ // 12C
		if(Ex>Ex_12C_s12) return 1; // s1/2-hole
		else return 2; //p3/2-hole
	}else if(Zt==8&&Nt==8){
		if(Ex>Ex_16O_s12) return 1;
		else if(Ex>Ex_16O_p32) return 2;
		else return 3;
	}else abort();
	return -1;
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
	if(nuc_daughter!=NULL) name_daughter = nuc_daughter->name;
	else{
		os.str("");
		os << Z_daughter+N_daughter << _nucleus_table->nuc_name[Z_daughter];
		name_daughter = os.str();
	}

	if(verbose>0){ 
		cout << "DecayMode(): Random = " << random << " : " << name_target.c_str() << " --> " << particle_name[decay_mode] << " + ";
		cout << name_daughter.c_str() << endl;
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
void Deexcitation::Decay(const bool breakflag)
/////////////////////////////////////////////
{
	cout << "Deexcitation::Decay()" << endl;
	
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
		cout << "mass_target   = " << mass_target << endl;
		cout << "mass_particle = " << mass_particle << endl;
		cout << "mass_daughter = " << mass_daughter << endl;
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
											particle_name[decay_mode],
											1,0,verbose);// trace flag==1, Ex_daughter==0
	double totalE_particle_bef = p_particle.totalE();
	p_particle.Boost(totalE_target,mom_target);// BOOST!
	double totalE_particle_aft = p_particle.totalE();

	Particle p_daughter(PDGion(Z_daughter,N_daughter),
											mass_daughter, // w/o excitation E
											mom_daughter,
											name_daughter,
											0,Ex_daughter,verbose); // intermediate state in default

	double totalE_daughter_bef = p_daughter.totalE();
	p_daughter.Boost(totalE_target,mom_target);
	double totalE_daughter_aft = p_daughter.totalE();

	// 
	double kE_target = totalE_target - mass_target;
	double totalE_ex_target = totalE_target + Ex_target; // w/ excitation E

	double totalE_bef = totalE_particle_bef + totalE_daughter_bef;
	double totalE_aft = totalE_particle_aft + totalE_daughter_aft;
	double totalE_ex_bef = totalE_bef + Ex_daughter; // w/ excitation E
	double totalE_ex_aft = totalE_aft + Ex_daughter;

	if(verbose>0){
		cout << "totalE_target = " << totalE_target << endl;
		cout << "kE_target = " << kE_target << endl;
		cout << "Ex_target = " << Ex_target << endl;
		cout << "totalE_ex_target = " << totalE_ex_target << endl;
		cout << "totalE_bef = " << totalE_bef << endl;
		cout << "totalE_aft = " << totalE_aft << endl;
		cout << "  diff = " << totalE_aft-totalE_bef << endl;
		cout << "totalE_ex_bef = " << totalE_ex_bef << endl;
		cout << "totalE_ex_aft = " << totalE_ex_aft << endl;
	}


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

	// Then, push back

	_particle->push_back(p_particle);
	// DoDeex loop will be end -> turn on track flag, because it is not intermediate state
	if(breakflag) p_daughter._flag=1;
	_particle->push_back(p_daughter);

	if(verbose>0){
		cout << "flag_particle = " << p_particle._flag << endl;
		cout << "flag_daughter = " << p_daughter._flag << endl;
	}
}

/////////////////////////////////////////////
double Deexcitation::ElementMassInMeV(TGeoElementRN* ele)
/////////////////////////////////////////////
{
	if(ele==0) return -1; // no profile can be found.
	double mass = ele->MassNo()*amu_c2
									 + ele->MassEx(); // (MeV)
	double mass_amu = mass/amu_c2;
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
bool Deexcitation::OpenROOT(const int Zt,const int Nt, const int Z, const int N, 
														const bool tree)
/////////////////////////////////////////////
{
	os.str("");
	os << getenv("TALYS_WORK") << "/output/";
	// single nucleon hole
	if(Zt+Nt==Z+N+1){
		if(Zt==6&&Nt==6) os << "12C/";
		else if(Zt==8&&Nt==8) os << "16O/";
		else return 0; // not supported
	}else{
		// multi-nucleon hole
		;
	}
	os << "Br_" << name_target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << ".root"; 
	//
	rootf = new TFile(os.str().c_str(),"READ");
	if(! rootf->IsOpen()) return 0;
	if(verbose>1){
		cout << "OpenRoot: " << os.str().c_str() << endl;
	}
	if(!tree) return 1;
	else return GetTTree(Z,N);
}

/////////////////////////////////////////////
bool Deexcitation::GetTTree(const int Z, const int N)
/////////////////////////////////////////////
{
	tree = (TTree*)rootf->Get("tree");
	if(tree==0) return 0; // not ttree
	tree->SetBranchAddress("Z",&_Z);
	tree->SetBranchAddress("N",&_N);
	tree->SetBranchAddress("Ex_bin",&_Ex_bin);
	tree->SetBranchAddress("Ex",&_Ex);
	tree->SetBranchAddress("Br",&_Br);
	tree->SetBranchAddress("REx_bin",&_REx_bin);
	tree->SetBranchAddress("REx",&_REx);
	tree->SetBranchAddress("RBr",&_RBr);
	return 1;
}


/////////////////////////////////////////////
bool Deexcitation::CreateTGraph(const int Z, const int N)
/////////////////////////////////////////////
{
	bool found=0;
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
		if(Z==_Z && N==_N){
			found=1;
			break;
		}
	}
	if(!found) return 0;

	for(int p=0;p<num_particle;p++){
		g_br[p] = new TGraph(_Ex_bin[p],_Ex[p],_Br[p]);
	}
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
		if(g_br[p]==0) return 0; // no tgraph
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
	if(point==g_br[mode]->GetN()) point--;
	if(point<0) point=0;
	if(verbose>0){
		cout << "GetBrExTGraph(): nearest_point = " << point << ",  diff_Ex = " << abs(ex-ex_t) << ", diff_ex(previous) = " << diff_ex << endl;
	}
	os.str("");
	os << "g_" << st.c_str() << "_br_ex_" << mode << "_" << point;
	g_br_ex = (TGraph*) rootf->Get(os.str().c_str());
	if(g_br_ex==0) return -1; // no tgraph

	return point;
}


/////////////////////////////////////////////
void Deexcitation::DeleteTGraphs()
/////////////////////////////////////////////
{
	for(int p=0;p<num_particle;p++){
		delete g_br[p];
		g_br[p]=0;
	}
	delete g_br_ex;
	g_br_ex=0;
}


/////////////////////////////////////////////
void Deexcitation::InitParticleVector()
/////////////////////////////////////////////
{
	if(_particle!=0){
		_particle->clear();
		delete _particle;
	}
	_particle = new vector<Particle>;
}
