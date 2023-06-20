#ifndef __NucleusTable__HH__
#define __NucleusTable__HH__

#include <map>
#include "Nucleus.hh"

using namespace std;

class NucleusTable{
  public:
	NucleusTable(){num_of_nuc=-1;};
  virtual ~NucleusTable(){;};

  bool ReadTables();
	int getID(const char* name);
	int GetNumofNuc(){return num_of_nuc;};
	
	Nucleus* GetNucleusPtr(const char* name);
	Nucleus* GetNucleusPtr(int id);
	Nucleus* GetNucleusPtr(int Z,int N);
	const char* nuc_name[9]
		= {"","H","He","Li","Be","B","C","N","O"}; // [Z]

  private:
	int num_of_nuc;
	Nucleus* _nucleus;
	map<string, int> _nucleus_id;
	map<string, int> :: iterator _p_id;
};
#endif
