#ifndef __NUCLEUS__HH__
#define __NUCLEUS__HH__

#include "consts.hh"

using namespace std;

class Nucleus{
  public :
	Nucleus();
  Nucleus(const char* Name,
          int z,
          int n);
	~Nucleus();

  char name[5];
  int Z;
  int N;
  int A;
  int id;

	bool flag_s; // 1 -> have separation energy file
	bool flag_target; // 1 -> target nucleus
	bool flag_data; // 1 -> have data in talys output (population)

	float sum_pop; // will be sum(total_pop[])

  float* total_pop;
  float** pop;
  float** Ex;
  int* Ex_bin;
	// [parity] [mother E bin]
	//		[parity] : 0 (negative), 1 (positive)

	// NOTATION
	// pop[p][mb][db]: population
	// Ex[p][mb][db]:  Daughter excitation energy
	// Ex_bin[p][mb]: number of daughter excitation energy bin
	// S[p]: Separation energy
	//		[particle][mother E bin][daughter bin]

	float***  pop_p;
	float*** Ex_p;
	int** Ex_bin_p;

	float* S; // separtion energy [particle]
	float min_S(); // minimum separation energy

	float GetPopP(int p,int mb); // return  sum_(daughter_bin) pop[p][mb][daughter_bin]

	bool CheckPop();
	

	private:
	void Init();
};
#endif
