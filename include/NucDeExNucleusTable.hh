#ifndef __NUCDEEXNUCLEUSTABLE__HH__
#define __NUCDEEXNUCLEUSTABLE__HH__

#include <map>

#include "NucDeExNucleus.hh"

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

class NucDeExNucleusTable{
  public:
  NucDeExNucleusTable();
#ifdef INCL_DEEXCITATION_NUCDEEX
  NucDeExNucleusTable(G4INCL::Config *config);
#endif
  virtual ~NucDeExNucleusTable(){;};

  bool ReadTables(const bool init_flag=1);
  int getID(const char* name);
  int GetNumofNuc(){return num_of_nuc;};
  void SetVerbose(const int v){verbose=v;};
  
  NucDeExNucleus* GetNucleusPtr(const char* name);
  NucDeExNucleus* GetNucleusPtr(int id);
  NucDeExNucleus* GetNucleusPtr(int Z,int N);
  NucDeExNucleus* GetNucleusPtrPDG(int PDG);
  const char* nuc_name[10]
    = {"","H","He","Li","Be","B","C","N","O","F"}; // [Z]

  private:
  int num_of_nuc;
  NucDeExNucleus* _nucleus;
  std::map<std::string, int> _nucleus_id;
  std::map<std::string, int> :: iterator _p_id;
  std::string PATH_NucDeEx_table;
  int verbose;
};
#endif
