#ifndef __NucleusTable__HH__
#define __NucleusTable__HH__

#include <map>
#include <Nucleus.hh>

using namespace std;

class NucleusTable{
  public:
	NucleusTable(){num_of_nuc=-1;};
  virtual ~NucleusTable(){;};

  bool ReadTables();
	int getID(const char* name);
	
	Nucleus* GetNucleusPtr(const char* name);
	Nucleus* GetNucleusPtr(int id);

  private:
	int num_of_nuc;
	Nucleus* _nucleus;
	map<string, int> _nucleus_id;
	map<string, int> :: iterator _p_id;

/*
  int getZ(const char* name);
  int getN(const char* name);
  int getA(const char* name);

  double getTotalPop(const char* name);
  void setTotalPop(const char* name, double pop);

  int getIndex(const char* name);
  void setIndex(const char* name, int i);

  const char* getName(int z, int n);

  int getExBin(const char* name);
  void setExBin(const char* name, int i);
  double getPop(const char* name, int i);
  void setPop(const char* name, int i, double pop);
  double getEx(const char* name, int i);
  void setEx(const char* name, int i, double Ex);
  
////////// Decay
  // gamma 
  int getExBinG(const char* name, int i);
  void setExBinG(const char*name, int i, int j);
  double getTotalPopG(const char* name, int i);
  const double ** getPopG(const char* name);
  double getPopG(const char* name, int i, int j);
  void setPopG(const char* name, int i, int j, double pop);
  const double ** getExG(const char* name);
  double getExG(const char* name, int i, int j);
  void setExG(const char* name, int i, int j, double Ex);
  
  // neutron
  int getExBinN(const char* name, int i);
  void setExBinN(const char*name, int i, int j);
  double getTotalPopN(const char* name, int i);
  const double ** getPopN(const char* name);
  double getPopN(const char* name, int i, int j);
  void setPopN(const char* name, int i, int j, double pop);
  const double ** getExN(const char* name);
  double getExN(const char* name, int i, int j);
  void setExN(const char* name, int i, int j, double Ex);

  // proton
  int getExBinP(const char* name, int i);
  void setExBinP(const char*name, int i, int j);
  double getTotalPopP(const char* name, int i);
  const double ** getPopP(const char* name);
  double getPopP(const char* name, int i, int j);
  void setPopP(const char* name, int i, int j, double pop);
  const double ** getExP(const char* name);
  double getExP(const char* name, int i, int j);
  void setExP(const char* name, int i, int j, double Ex);

  // alpha
  int getExBinA(const char* name, int i);
  void setExBinA(const char*name, int i, int j);
  double getTotalPopA(const char* name, int i);
  const double ** getPopA(const char* name);
  double getPopA(const char* name, int i, int j);
  void setPopA(const char* name, int i, int j, double pop);
  const double ** getExA(const char* name);
  double getExA(const char* name, int i, int j);
  void setExA(const char* name, int i, int j, double Ex);

  // deuteron
  int getExBinD(const char* name, int i);
  void setExBinD(const char*name, int i, int j);
  double getTotalPopD(const char* name, int i);
  const double ** getPopD(const char* name);
  double getPopD(const char* name, int i, int j);
  void setPopD(const char* name, int i, int j, double pop);
  const double ** getExD(const char* name);
  double getExD(const char* name, int i, int j);
  void setExD(const char* name, int i, int j, double Ex);

  // triton
  int getExBinT(const char* name, int i);
  void setExBinT(const char*name, int i, int j);
  double getTotalPopT(const char* name, int i);
  const double ** getPopT(const char* name);
  double getPopT(const char* name, int i, int j);
  void setPopT(const char* name, int i, int j, double pop);
  const double ** getExT(const char* name);
  double getExT(const char* name, int i, int j);
  void setExT(const char* name, int i, int j, double Ex);

  // he3
  int getExBinH(const char* name, int i);
  void setExBinH(const char*name, int i, int j);
  double getTotalPopH(const char* name, int i);
  const double ** getPopH(const char* name);
  double getPopH(const char* name, int i, int j);
  void setPopH(const char* name, int i, int j, double pop);
  const double ** getExH(const char* name);
  double getExH(const char* name, int i, int j);
  void setExH(const char* name, int i, int j, double Ex);
*/

/*
  map<string, nucleus> _nucleus_table;
  map<string, nucleus> :: iterator _p;
*/
};

#endif
