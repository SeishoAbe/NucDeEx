#ifndef __NUCLEUS__HH__
#define __NUCLEUS__HH__

using namespace std;

class nucleus{
  public :
  nucleus(const char* name,
          int z,
          int n);
  virtual ~nucleus(){}

  char Name[5];
  int Z;
  int N;
  int A;
  double total_pop;
  double* pop;
  double* Ex;
  int Ex_bin;
  int index;

	// NOTATION
	// pop_x[][]: population
	// Ex_x[][]:  Excitation energy
	// Ex_bin_x[]: Excitatioln energy bin
	// S_x: Separation energy
  
  // gamma
  double** pop_g;
  double** Ex_g;
  int* Ex_bin_g;
	double S_g;

  // neutron
  double** pop_n;
  double** Ex_n;
  int* Ex_bin_n;
	double S_n;

  // proton
  double** pop_p;
  double** Ex_p;
  int* Ex_bin_p;
	double S_p;

  // alpha
  double** pop_a;
  double** Ex_a;
  int* Ex_bin_a;
	double S_a;

  // deuteron
  double** pop_d;
  double** Ex_d;
  int* Ex_bin_d;
	double S_d;

  // triton
  double** pop_t;
  double** Ex_t;
  int* Ex_bin_t;
	double S_t;

  // he3
  double** pop_h;
  double** Ex_h;
  int* Ex_bin_h;
	double S_h;

  const int array=100;
};
#endif
