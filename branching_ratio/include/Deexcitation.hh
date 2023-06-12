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
	// Simulation method called by DoDeex()
	int DecayMode(double Ex);
	int NearestExPoint(double Ex, int decay_mode);
	bool DaughterExPoint(double &d_Ex, int &d_point);

	// ROOT related methods & members
	bool ReadROOT(const char* name);
	TFile* rootf;
	TTree* tree;
	TGraph* g_br[num_particle];
	TGraph* g_br_ex;

	// Target info
	int Z_target, N_target;
	double Ex_target;
	string name_target;

	// nucleus table
	NucleusTable* _nucleus_table;

	// Random generator
	TRandom3* rndm;

	// others
	int verbose;
	ostringstream os;
};
#endif
