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

	bool flag_s;
	// 0 -> does not have separation energy file
	// 1 -> HAVE IT

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
	

	private:
	void Init();
};
#endif
