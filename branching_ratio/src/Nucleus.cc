#include "Nucleus.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>

using namespace std;


nucleus::nucleus(const char* name, int z, int n){
  strcpy(Name,name);
  Z=z;
  N=n;
  A=z+n;
  total_pop=0.0;

  pop = new double[array];
  Ex = new double[array];
  Ex_bin=0.0;

  pop_g = new double*[array];
  pop_n = new double*[array];
  pop_p = new double*[array];
  pop_a = new double*[array];
  pop_d = new double*[array];
  pop_t = new double*[array];
  pop_h = new double*[array];
  Ex_g = new double*[array];
  Ex_n = new double*[array];
  Ex_p = new double*[array];
  Ex_a = new double*[array];
  Ex_d = new double*[array];
  Ex_t = new double*[array];
  Ex_h = new double*[array];
  Ex_bin_g = new int[array];
  Ex_bin_n = new int[array];
  Ex_bin_p = new int[array];
  Ex_bin_a = new int[array];
  Ex_bin_d = new int[array];
  Ex_bin_t = new int[array];
  Ex_bin_h = new int[array];
  for(int i=0;i<array;i++){
    pop[i]=0.;
    Ex[i]=0.;
    pop_g[i] = new double[array];
    pop_n[i] = new double[array];
    pop_p[i] = new double[array];
    pop_a[i] = new double[array];
    pop_d[i] = new double[array];
    pop_t[i] = new double[array];
    pop_h[i] = new double[array];
    Ex_g[i] = new double[array];
    Ex_n[i] = new double[array];
    Ex_p[i] = new double[array];
    Ex_a[i] = new double[array];
    Ex_d[i] = new double[array];
    Ex_t[i] = new double[array];
    Ex_h[i] = new double[array];
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
  index=0;
	S_g=S_n=S_p=S_a=S_d=S_t=S_h=0.;
}



