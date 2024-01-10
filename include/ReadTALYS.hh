#ifndef __READTALYS__HH__
#define __READTALYS__HH__

#include <string>
#include <ostream>

#include "NucDeExNucleusTable.hh"
#include "NucDeExConsts.hh"

class ReadTALYS{
  public: 
  //ReadTALYS(){;};
  ReadTALYS(const char* filename, NucDeExNucleusTable* nuc);
  virtual ~ReadTALYS(){;};

  bool Read();
  void SetVerboseLevel(int v){ _verbose = v;};

  private:
  std::string _filename;
  NucDeExNucleusTable* _nucleus_table;
  std::ifstream* _ifs;
  int _verbose;

  const float check_criteria=0.05;

  void SetKeywords();

  std::string* keyword_population;
  std::string* keyword_N;
  std::string* keyword_parity;
  std::string* keyword_before_decay;
  std::string* keyword_decay;
  std::string* keyword_total;
  std::string* keyword_bin_mother;
  std::string* keyword_parity_mother;
  std::string* keyword_parity_daughter;
  std::string* keyword_discrete;
  std::string* keyword_discrete_br;

  std::ostringstream* os;
};

#endif
