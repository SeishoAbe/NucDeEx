#ifndef __READTALYS__HH__
#define __READTALYS__HH__

#include <string>
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

	string keyword_multiple_emission;
	string keyword_population;
	string keyword_N;
	string keyword_before_decay;
};

#endif
