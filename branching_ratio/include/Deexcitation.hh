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
	virtual ~Deexcitation(){;};

	void DoDeex(int Z,int N, double Ex);

	void SetSeed(int s){ ran->SetSeed(s) ;};
	int GetSeed(){ran->GetSeed();};
	void SetVerbose(int v){ verbose=v; };

	private:
	// root 
	bool ReadROOT();
	TFile* rootf;
	TTree* tree;
	TGraph* g_br[num_particle];
	TGraph* g_br_ex;

	// target info
	int Z_target, N_target;
	string name_target;

	NucleusTable* _nucleus_table;

	// rand
	TRandom3* ran;

	// others
	int verbose;
	ostringstream os;
};
#endif
