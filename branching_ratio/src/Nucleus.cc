#include "Nucleus.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include "consts.hh"

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
	id=0;
	sum_pop = 0;
	total_pop = new float[parity];
  pop = new float*[parity];
  Ex = new float*[parity];
	Ex_bin = new int[parity];
	for(int i=0;i<parity;i++){
		total_pop[i] = 0;
		Ex_bin[i] = 0;
		pop[i] = new float[bins];
		Ex[i] = new float[bins];
		for(int j=0;j<bins;j++){
			Ex[i][j] = -1;
		}
	}

	pop_p = new float**[num_particle];
	Ex_p = new float**[num_particle];
	Ex_bin_p = new int*[num_particle];
	S = new float[num_particle];

	for(int i=0;i<num_particle;i++){
		S[i]=0;
		pop_p[i] = new float*[bins];
		Ex_p[i] = new float*[bins];
		Ex_bin_p[i] = new int[bins];
		for(int j=0;j<bins;j++){
			pop_p[i][j] = new float[bins];
			Ex_p[i][j] = new float[bins];
		}

	}
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

  for(int i=0;i<bins;i++){
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
