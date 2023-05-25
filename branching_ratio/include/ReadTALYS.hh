#ifndef __READTALYS__HH__
#define __READTALYS__HH__

#include <string>
#include <ostream>
#include "NucleusTable.hh"
#include "consts.hh"

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

	string* keyword_multiple_emission;
	string* keyword_population;
	string* keyword_N;
	string* keyword_parity;
	string* keyword_before_decay;
	string* keyword_decay;
	string* keyword_total;
	string* keyword_bin_mother;
	string* keyword_parity_mother;
	string* keyword_parity_daughter;

	ostringstream* os;
};

#endif
