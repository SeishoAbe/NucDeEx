#ifndef __NUCLEUS__HH__
#define __NUCLEUS__HH__

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
  float total_pop;
  float* pop;
  float* Ex;
  int Ex_bin;
  int id;

	// NOTATION
	// pop_x[][]: population
	// Ex_x[][]:  Excitation energy
	// Ex_bin_x[]: Excitatioln energy bin
	// Sx: Separation energy
  
  // gamma
  float** pop_g;
  float** Ex_g;
  int* Ex_bin_g;
	float Sg;

  // neutron
  float** pop_n;
  float** Ex_n;
  int* Ex_bin_n;
	float Sn;

  // proton
  float** pop_p;
  float** Ex_p;
  int* Ex_bin_p;
	float Sp;

  // alpha
  float** pop_a;
  float** Ex_a;
  int* Ex_bin_a;
	float Sa;

  // deuteron
  float** pop_d;
  float** Ex_d;
  int* Ex_bin_d;
	float Sd;

  // triton
  float** pop_t;
  float** Ex_t;
  int* Ex_bin_t;
	float St;

  // he3
  float** pop_h;
  float** Ex_h;
  int* Ex_bin_h;
	float Sh;
	
  const int array=100;
	bool flag_s;
	// 0 -> does not have sep E file
	// 1 -> HAVE IT

	private:
	void Init();
};
#endif
