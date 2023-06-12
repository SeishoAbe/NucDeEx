#ifndef __DEEXCITATION__HH__
#define __DEEXCITATION__HH__

#include <string>
#include <ostream>
#include "NucleusTable.hh"
#include "consts.hh"

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

	void DoDeex(int Z,int N, double Ex);

	void SetSeed(int s){ rndm->SetSeed(s) ;};
	int  GetSeed(){return rndm->GetSeed();};
	void SetVerbose(int v){ verbose=v; };

	private:
	// --- Simulation method called by DoDeex()
	int DecayMode(double Ex);
	int NearestExPoint(double Ex, int decay_mode);
	bool DaughterExPoint(double &d_Ex, int &d_point);
	vector<TParticle*> Decay();
	const char* PDGion(int Z,int N);

	// --- ROOT related methods & members
	bool OpenROOT(const char* name);
	bool GetBrTGraph();
	TFile* rootf;
	TTree* tree;
	TGraph* g_br[num_particle];
	TGraph* g_br_ex;
	TRandom3* rndm;
	TDatabasePDG* pdg;
	TGeoManager* geo;
	TGeoElementTable* element_table;

	double ElementMassInMeV(TGeoElementRN* ele);

	// -- Nucleus info
	// target info
	int Z_target, N_target;
	double Ex_target;
	Nucleus* nuc_target;
	double mass_target;
	string name_target;

	// daughter info
	int Z_daughter, N_daughter;
	double Ex_daughter;
	Nucleus* nuc_daughter;
	double mass_daughter;

	// particle info
	double mass_particle;
	double S;
	double Qvalue;


	// others
	NucleusTable* _nucleus_table;
	int verbose;
	ostringstream os;
};
#endif
