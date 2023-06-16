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
	Deexcitation();
	virtual ~Deexcitation();

	void DoDeex(const int Z, const int N, const double Ex,
							const TVector3& mom=TVector3(0,0,0));

	void SetSeed(int s){ rndm->SetSeed(s) ;};
	int  GetSeed(){return rndm->GetSeed();};
	TRandom3* GetTRandom3(){ return rndm; };
	void SetVerbose(int v){ verbose=v; };
	void SetEventID(int id){ eventID=id;};
	int  GetEventID(){ return eventID; };
	vector<Particle> GetParticleVector(){return _particle;};

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
	TTree* tree;
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
	//double kE_target;
	//int PDG_target;
	Nucleus* nuc_target;
	string name_target;

	// daughter nucleus info
	int Z_daughter, N_daughter;
	double Ex_daughter;
	double mass_daughter;
	TVector3 mom_daughter;
	//double kE_daughter;
	//int PDG_daughter;
	Nucleus* nuc_daughter;
	string name_daughter;

	// decay particle info
	double mass_particle;
	TVector3 mom_particle;
	//double kE_particle;

	// for output
	void InitParticleVector();
	vector<Particle> _particle;
	int PDGion(int Z,int N);

	// others
	NucleusTable* _nucleus_table;
	int verbose;
	int eventID;
	ostringstream os;
	const double check_criteria=1e-3;
};
#endif
