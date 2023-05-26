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
	flag_s=flag_target=flag_data=0;
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
			pop[i][j] = 0;
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
			Ex_bin_p[i][j]=0;
			pop_p[i][j] = new float[bins];
			Ex_p[i][j] = new float[bins];
			for(int k=0;k<bins;k++){
				pop_p[i][j][k] = Ex_p[i][j][k] =0;
			}
		}

	}
}
//////////////////
float Nucleus::min_S()
//////////////////
{
	float min=1e9;
	for(int i=0;i<num_particle;i++){
		if(i==0) continue; // remove gamma
		if(min>S[i]) min=S[i];
	}
	return min;
}


//////////////////
float Nucleus::GetPopDaughterBinSum(int p,int mb)
//////////////////
{
	float population=0;
	//for(int i=0;i<=Ex_bin_p[p][mb];i++){
	for(int i=0;i<bins;i++){
		population += pop_p[p][mb][i];
	}

	return population;
}

//////////////////
float Nucleus::GetPopParticleDaughterBinSum(int mb)
//////////////////
{
	float population=0;
	for(int p=0;p<num_particle;p++){
		population += GetPopDaughterBinSum(p,mb);
	}
	return population;
}

//////////////////
bool Nucleus::CheckTotalPop()
//////////////////
{
	float sum_pop_check=0;
	for(int par=0;par<parity;par++){// parity loop
		float total_pop_check=0;
		//for(int i=0;i<=Ex_bin[par];i++){ // ex bin loop
		for(int i=0;i<bins;i++){
			total_pop_check += pop[par][i];
		}
		if(!(total_pop[par]>0)) continue;
		if( abs(total_pop[par]-total_pop_check)/total_pop[par] > 0.01 ){
			cerr << "ERROR: pality dependent total population is not reproduced" << endl;
			cerr << name << " parity=" << par << endl;
			cerr << "total population= " << total_pop[par] << "  summed population=" << total_pop_check << endl;
			return 0;
		}
		sum_pop_check += total_pop[par];
	}

	if( abs(sum_pop-sum_pop_check)/sum_pop > 0.01 ){
		cerr << "ERROR: pality sum total population is not reproduced" << endl;
		cerr << name << endl;
		cerr << "sum population= " << sum_pop << "  summed sum population=" << sum_pop_check << endl;
		return 0;
	}
}

//////////////////
bool Nucleus::CheckPop()
//////////////////
{
	for(int i=0;i<Ex_bin[0];i++){ // ex bin loop
		if(Ex[0][i]<min_S()) continue;
		float population=0;
		for(int par=0;par<parity;par++){ // parity loop
			population += pop[par][i]; // population obtained from "population"
		}
		if(!(population>0)) continue; // skip

		// sum population for daughter 
		float population_check=GetPopParticleDaughterBinSum(i);
		
		if( abs(population_check-population)/population > 0.01 ){
			cerr << "ERROR: population is not reproduced" << endl;
			cerr << name << " bin=" << i << endl;
			cerr << "population= " << population << "  summed population=" << population_check << endl;
			return 0;
		}
	}

	cout << "Population for " << name << " looks OK" << endl;
	return 1;
}

//////////////////
bool Nucleus::CheckEx()
//////////////////
{
	if(Ex_bin[0]!=Ex_bin[1]){
		cerr << "ERROR: Ex_bin[parity] are different" << endl;
		cerr << Ex_bin[0] << " " << Ex_bin[1] << endl;
		return 0;
	}

	for(int i=0;i<bins;i++){
		if(Ex[0][i] != Ex[1][i]){
			return 0;
		}
	}
	return 1;
}

//////////////////
Nucleus::~Nucleus(){
//////////////////
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
