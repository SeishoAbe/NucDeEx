#include "Nucleus.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>

using namespace std;

///////////////
Nucleus::Nucleus()
///////////////
{
	Init();
}

///////////////
Nucleus::Nucleus(const char* Name, int z, int n)
///////////////
{
  strcpy(name,Name);
  Z=z;
  N=n;
  A=z+n;
	Init();
}

///////////////
void Nucleus::Init()
///////////////
{
	flag_s=0;
	total_pop=0.0;

  pop = new float[array];
  Ex = new float[array];
  Ex_bin=0;

  pop_g = new float*[array];
  pop_n = new float*[array];
  pop_p = new float*[array];
  pop_a = new float*[array];
  pop_d = new float*[array];
  pop_t = new float*[array];
  pop_h = new float*[array];
  Ex_g = new float*[array];
  Ex_n = new float*[array];
  Ex_p = new float*[array];
  Ex_a = new float*[array];
  Ex_d = new float*[array];
  Ex_t = new float*[array];
  Ex_h = new float*[array];
  Ex_bin_g = new int[array];
  Ex_bin_n = new int[array];
  Ex_bin_p = new int[array];
  Ex_bin_a = new int[array];
  Ex_bin_d = new int[array];
  Ex_bin_t = new int[array];
  Ex_bin_h = new int[array];
  for(int i=0;i<array;i++){
    pop[i]=0;
    Ex[i]=-1;
    pop_g[i] = new float[array];
    pop_n[i] = new float[array];
    pop_p[i] = new float[array];
    pop_a[i] = new float[array];
    pop_d[i] = new float[array];
    pop_t[i] = new float[array];
    pop_h[i] = new float[array];
    Ex_g[i] = new float[array];
    Ex_n[i] = new float[array];
    Ex_p[i] = new float[array];
    Ex_a[i] = new float[array];
    Ex_d[i] = new float[array];
    Ex_t[i] = new float[array];
    Ex_h[i] = new float[array];
    Ex_bin_g[i] = 0;
    Ex_bin_n[i] = 0;
    Ex_bin_p[i] = 0;
    Ex_bin_a[i] = 0;
    Ex_bin_d[i] = 0;
    Ex_bin_t[i] = 0;
    Ex_bin_h[i] = 0;
    for(int j=0;j<array;j++){
      pop_g[i][j] = pop_n[i][j] = pop_p[i][j] = pop_a[i][j]
        = pop_d[i][j] = pop_t[i][j] = pop_h[i][j] = 0.;
      Ex_g[i][j] = Ex_n[i][j] = Ex_p[i][j] = Ex_a[i][j] 
        = Ex_d[i][j] = Ex_t[i][j] = Ex_h[i][j] = 0.;
    }
  }
  id=0;
	Sg=Sn=Sp=Sa=Sd=St=Sh=-1;
}


Nucleus::~Nucleus(){
	/*
	delete[] pop;
	delete[] Ex;

	delete[] Ex_bin_g;
	delete[] Ex_bin_n;
	delete[] Ex_bin_p;
	delete[] Ex_bin_a;
	delete[] Ex_bin_d;
	delete[] Ex_bin_t;
	delete[] Ex_bin_h;

  for(int i=0;i<array;i++){
		delete[] pop_g[i];
		delete[] pop_n[i];
		delete[] pop_p[i];
		delete[] pop_a[i];
		delete[] pop_d[i];
		delete[] pop_t[i];
		delete[] pop_h[i];
		delete[] Ex_g[i];
		delete[] Ex_n[i];
		delete[] Ex_p[i];
		delete[] Ex_a[i];
		delete[] Ex_d[i];
		delete[] Ex_t[i];
		delete[] Ex_h[i];
	}

	delete[] pop_g;
	delete[] pop_n;
	delete[] pop_p;
	delete[] pop_a;
	delete[] pop_d;
	delete[] pop_t;
	delete[] pop_h;
	delete[] Ex_g;
	delete[] Ex_n;
	delete[] Ex_p;
	delete[] Ex_a;
	delete[] Ex_d;
	delete[] Ex_t;
	delete[] Ex_h;
	*/
	;
}
