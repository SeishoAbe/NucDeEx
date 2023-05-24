#ifndef __NUCLEUS__HH__
#define __NUCLEUS__HH__

using namespace std;

class Nucleus{
  public :
	Nucleus();
  Nucleus(const char* name,
          int z,
          int n);
	~Nucleus();

  char Name[5];
  int Z;
  int N;
  int A;
  double total_pop;
  double* pop;
  double* Ex;
  int Ex_bin;
  int id;

	// NOTATION
	// pop_x[][]: population
	// Ex_x[][]:  Excitation energy
	// Ex_bin_x[]: Excitatioln energy bin
	// Sx: Separation energy
  
  // gamma
  double** pop_g;
  double** Ex_g;
  int* Ex_bin_g;
	double Sg;

  // neutron
  double** pop_n;
  double** Ex_n;
  int* Ex_bin_n;
	double Sn;

  // proton
  double** pop_p;
  double** Ex_p;
  int* Ex_bin_p;
	double Sp;

  // alpha
  double** pop_a;
  double** Ex_a;
  int* Ex_bin_a;
	double Sa;

  // deuteron
  double** pop_d;
  double** Ex_d;
  int* Ex_bin_d;
	double Sd;

  // triton
  double** pop_t;
  double** Ex_t;
  int* Ex_bin_t;
	double St;

  // he3
  double** pop_h;
  double** Ex_h;
  int* Ex_bin_h;
	double Sh;
	
  const int array=100;
	bool flag_s;
	// 0 -> does not have sep E file
	// 1 -> HAVE IT

	private:
	void Init();
};
#endif
