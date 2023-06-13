#ifndef __DEEXCITATION__HH__
#define __DEEXCITATION__HH__

#include <string>
#include <ostream>
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
							const TVector3* dir=0);

	void SetSeed(int s){ rndm->SetSeed(s) ;};
	int  GetSeed(){return rndm->GetSeed();};
	void SetVerbose(int v){ verbose=v; };

	private:
	// --- Simulation method called by DoDeex() --- //
	int DecayMode(double Ex);
	bool DaughterExPoint(double *d_Ex, int *d_point);
	void Decay();


	// --- ROOT related methods & members --- //
	bool OpenROOT(const char* name);
	bool GetBrTGraph(string st);
	int  GetBrExTGraph(string st, double ex_t, int mode); 
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
	// target nucleus info
	int Z_target, N_target;
	double Ex_target;
	Nucleus* nuc_target;
	double mass_target;
	string name_target;
	TVector3* dir_target;

	// daughter nucleus info
	int Z_daughter, N_daughter;
	double Ex_daughter;
	Nucleus* nuc_daughter;
	double mass_daughter;

	// decay particle info
	double mass_particle;
	double S;
	double Qvalue;

	// for output
	vector<Particle*> _particle;
	const char* PDGion(int Z,int N);

	// others
	NucleusTable* _nucleus_table;
	int verbose;
	ostringstream os;
};
#endif
