#ifndef __NucleusTable__HH__
#define __NucleusTable__HH__

#include <map>
#include "Nucleus.hh"

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

class NucleusTable{
  public:
  NucleusTable();
#ifdef INCL_DEEXCITATION_NUCDEEX
  NucleusTable(G4INCL::Config *config);
#endif
  virtual ~NucleusTable(){;};

  bool ReadTables(const bool init_flag=1);
	int getID(const char* name);
	int GetNumofNuc(){return num_of_nuc;};
  void SetVerbose(const int v){verbose=v;};
	
	Nucleus* GetNucleusPtr(const char* name);
	Nucleus* GetNucleusPtr(int id);
	Nucleus* GetNucleusPtr(int Z,int N);
	Nucleus* GetNucleusPtrPDG(int PDG);
	const char* nuc_name[10]
		= {"","H","He","Li","Be","B","C","N","O","F"}; // [Z]

  private:
	int num_of_nuc;
	Nucleus* _nucleus;
	std::map<std::string, int> _nucleus_id;
	std::map<std::string, int> :: iterator _p_id;
  std::string PATH_NucDeEx_table;
  int verbose;
};
#endif
