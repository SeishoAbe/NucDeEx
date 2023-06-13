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
	int maxlevelsbin;

	bool flag_s; // 1 -> have separation energy file
	bool flag_target; // 1 -> target nucleus
	bool flag_data; // 1 -> have data in talys output (population)

	float** level_br;
	float* level_Ex;
	int level_Ex_bin;
	// [Ex bin][daughter Ex bin]
	// For (gamma) level info

	float sum_pop; // will be sum(total_pop[])

  float* total_pop;
  float** pop;
  float** Ex;
  int* Ex_bin;
	// [parity] [mother E bin]
	//		[parity] : 0 (negative), 1 (positive)
	bool* flag_decay_data; // does it have decay data?
	 //[mother E bin]

	// NOTATION
	// pop[p][mb][db]: population
	// Ex[p][mb][db]:  Daughter excitation energy
	// Ex_bin[p][mb]: number of daughter excitation energy bin
	//		[particle][mother E bin][daughter bin]

	float***  pop_p;
	float*** Ex_p;
	int** Ex_bin_p;

	float* S; // separtion energy [particle]
	float min_S(); // minimum separation energy

	float GetPopDaughterBinSum(int p,int mb); // return  sum_(daughter_bin) pop[p][mb][daughter_bin]
		// (particle, mother ex bin)
	float GetPopParticleDaughterBinSum(int mb);

	float GetPopParitySum(int mb){ return pop[0][mb] + pop[1][mb];};

	bool CheckTotalPop();
	bool CheckPop();
	//bool CheckPop(int i);
	bool CheckEx();
	
	private:
	void Init();
	const float check_criteria=0.05;
};
#endif
