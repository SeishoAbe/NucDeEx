#ifndef __READTALYS__HH__
#define __READTALYS__HH__

#include <string>
#include <ostream>
#include "NucleusTable.hh"

using namespace std;

class ReadTALYS{
	public: 
	//ReadTALYS(){;};
	ReadTALYS(const char* filename,NucleusTable* nuc);
	virtual ~ReadTALYS(){;};

	bool Read();

	private:
	string _filename;
	NucleusTable* _nucleus_table;
	ifstream* _ifs;

	void SetKeywords();

	bool ReadPopulation();
	bool ReadDecay();

	string* keyword_multiple_emission;
	string* keyword_population;
	string* keyword_N;
	string* keyword_parity;
	string* keyword_before_decay;
	string* keyword_decay;
	string* keyword_bin_mother;
	int skip_after_population=4;

	const int num_particle=7;
	string*  particle_name;

	ostringstream* os;
};

#endif
