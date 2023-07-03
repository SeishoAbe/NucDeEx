#ifndef __DEEXCITATION__HH__
#define __DEEXCITATION__HH__

#include <string>
#include <ostream>
#include <vector>
#include "consts.hh"
#include "NucleusTable.hh"
#include "Particle.hh"

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TGeoManager.h>
#include <TGeoElement.h>

using namespace std;

class Deexcitation{
	public:
	Deexcitation(const int ldmodel=1, const bool parity_optmodall=1);
	virtual ~Deexcitation();

	// --- Main function ---//
	int DoDeex(const int Zt, const int Nt,
		 				 const int Z, const int N, const int shell, const double Ex,
						 const TVector3& mom=TVector3(0,0,0));
	// Zt, Nt  : Target nucleus Z and N (supports 16O or 12C currently)
	// Z, N    : Residual nucleus Z and N (having excitation energy)
	// shell  0 : Shell level will be determined according to Ex (box cut)
	//        1 : s1/2-hole
	//			  2 : p3/2-hole
	//        3 : p1/2-hole (only for 16O target)
	// Ex      : Only used if shell==1
	// mom     : 3D momentum of residual nucleus
	// 
	// return  0 : The target or residual nucleus is not supported
	//         1 : Success
	//        -1 : Fatal error

	// ---Sub functions --- //
	// supports multi-nucleon disapperance
	int DoDeex_talys(const int Zt, const int Nt,
							     const int Z, const int N, const double Ex,
							     const TVector3& mom=TVector3(0,0,0));
	// return 0: No talys tables. Nothing to do.
	//				1: Sucess
	//       -1: Very rarely happens due to production of 5H etc. (no mass table)
	//           If you find this, please execute this again

	// only for single-nucleon disappearance
	int DoDeex_p32(const int Zt, const int Nt,
								 const int Z, const int N,
							   const TVector3& mom=TVector3(0,0,0));
	// Note: Energy level is determined irrelevant to Ex
	// return  1: Sucess
	//        -1 : Fatal error
	
	void AddGSNucleus(const int Z, const int N, const TVector3& mom=TVector3(0,0,0));
	// Just add g.s. nucleus to the vector.
	// this function will be used in p1/2-hole, etc.
	int ExtoShell(const int Zt, const int Nt, const double Ex);

	void SetSeed(int s){ rndm->SetSeed(s) ;};
	int  GetSeed(){return rndm->GetSeed();};
	TRandom3* GetTRandom3(){ return rndm; };
	void SetVerbose(int v){ verbose=v; };
	void SetEventID(int id){ eventID=id;};
	int  GetEventID(){ return eventID; };
	vector<Particle>* GetParticleVector(){return _particle;};
	NucleusTable* GetNucleusTablePtr(){ return _nucleus_table;};
	int GetShell(){return _shell;};

	private:
	// --- Simulation method called by DoDeex() --- //
	int DecayMode(const double Ex);
	bool DaughterExPoint(double *d_Ex, int *d_point); //call by pointer
	void Decay(const bool breakflag);

	// --- ROOT related methods & members --- //
	bool OpenROOT(const char* name);
	bool GetBrTGraph(const string st);
	int  GetBrExTGraph(const string st, const double ex_t, const int mode); 
		// The nearest TGraph point will be returned
	TFile* rootf;
	//TTree* tree;
	TGraph* g_br[num_particle];
	TGraph* g_br_ex;
	TRandom3* rndm;
	TDatabasePDG* pdg;
	TGeoManager* geo;
	TGeoElementTable* element_table;
	double ElementMassInMeV(TGeoElementRN* ele);

	// --- Decay information (in MeV) --- //
	int decay_mode;
	double S;
	double Qvalue;

	// target nucleus info
	int Z_target, N_target;
	double Ex_target;
	double mass_target;
	TVector3 mom_target;
	Nucleus* nuc_target;
	string name_target;

	// daughter nucleus info
	int Z_daughter, N_daughter;
	double Ex_daughter;
	double mass_daughter;
	TVector3 mom_daughter;
	Nucleus* nuc_daughter;
	string name_daughter;

	// decay particle info
	double mass_particle;
	TVector3 mom_particle;

	// for output
	void InitParticleVector();
	vector<Particle> *_particle;
	int PDGion(int Z,int N);

	// others
	NucleusTable* _nucleus_table;
	int ldmodel;
	bool parity_optmodall;
	int verbose;
	int eventID;
	ostringstream os;
	const double check_criteria=1e-3;
	int _shell;

	// Constatnts for Ex to shell
	const double Ex_12C_s12=16.0;
	const double Ex_16O_s12=14.0;
	const double Ex_16O_p32=4.0;

	// Constants for (p3/2)-1 Br
	// 11B* (Panin et al., Phys. Lett. B 753 204-210. Experimental data)
	static const int Nlevel_p32_11B = 3; 
	const double E_p32_11B[Nlevel_p32_11B]={0.,2.125,5.020};
	const double Br_p32_11B[Nlevel_p32_11B]={0.82,0.10,0.08};
	// 11C* (assume analogy of 11B*)
	//      the same Br, but energy is different
	static const int Nlevel_p32_11C = 3;
	const double E_p32_11C[Nlevel_p32_11C]={0.,2.000,4.804};
	const double Br_p32_11C[Nlevel_p32_11C]={0.82,0.10,0.08};

	// 15N* (Ejili, Phys. Rev. C 58, 3)
	//		 Deexcitation from 9.93 MeV will be described from TALYS data
	//     Deexcitation from 10.70 is taken from Ejiri (100% proton)
	static const int Nlevel_p32_15N = 3;
	const double E_p32_15N[Nlevel_p32_15N]={6.32,9.93,10.70};
	const double Br_p32_15N[Nlevel_p32_15N]={0.872,0.064,0.064};
	// 15O* (Ejiri , Phys. Rev. C 58, 3)
	static const int Nlevel_p32_15O=3;
	const double E_p32_15O[Nlevel_p32_15O]={6.18,9.61,10.48};
	const double Br_p32_15O[Nlevel_p32_15O]={0.872,0.064,0.064}; // guess
	//const double Br_p32_15O[Nlevel_p32_15O]={1.,0,0}; // original Ejiri's value
};
#endif
