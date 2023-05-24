#include "NucleusTable.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>

using namespace std;

NucleusTable* NucleusTable::_NucleusTable = new NucleusTable();

///////////////

bool NucleusTable::ReadTable(const char* filename){
  std::ifstream ifs(filename);
  if(!ifs.is_open()){
    std::cerr << "ERROR : Cannot open " << filename << std::endl;
    abort;
  }
  char buf[500];
  //const char* name;
  char name[5];
  int z,n;
  while(ifs.getline(buf,sizeof(buf))){
    if(buf[0]=='#') continue;
    std::istringstream(buf) >> name >> z >> n;
    _nucleus_table.insert(make_pair(name,nucleus(name,z,n)));
  }
  ifs.close();
  return true;
}

int NucleusTable::getZ(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Z);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
int NucleusTable::getN(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).N);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
int NucleusTable::getA(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).A);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getTotalPop(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).total_pop);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setTotalPop(const char* name, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).total_pop = pop;
  }
}
int NucleusTable::getIndex(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).index);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setIndex(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).index = i;
  }
}
int NucleusTable::getExBin(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBin(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin = i;
  }
}
//
double NucleusTable::getPop(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPop(const char* name, int i, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop[i] = pop;
  }
}
double NucleusTable::getEx(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setEx(const char* name, int i, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex[i] = Ex;
  }
}


const char* NucleusTable::getName(int z, int n){
  for(auto i=_nucleus_table.begin(); i!=_nucleus_table.end() ; i++){
    if(int((i->second).Z) == z && int((i->second).N)==n){
      return (const char*) (i->second).Name;
    }
  }
  std::cerr << "ERROR : Cannnot find such nucleus Z " << z << " N " << n<< std::endl;
  abort();
}


////////// Decay
// gamma
int NucleusTable::getExBinG(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_g[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinG(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_g[i] = j;
  }
}
double NucleusTable::getTotalPopG(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_g =  (const double**)((_p->second).pop_g);
    double total_pop_g=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_g += pop_g[i][j];
    }
    return total_pop_g;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopG(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_g);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopG(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_g[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopG(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_g[i][j] = pop;
  }
}
const double** NucleusTable::getExG(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_g);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExG(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_g[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExG(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_g[i][j] = Ex;
  }
}

// neutron
int NucleusTable::getExBinN(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_n[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinN(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_n[i] = j;
  }
}
double NucleusTable::getTotalPopN(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_n =  (const double**)((_p->second).pop_n);
    double total_pop_n=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_n += pop_n[i][j];
    }
    return total_pop_n;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopN(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_n);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopN(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_n[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopN(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_n[i][j] = pop;
  }
}
const double** NucleusTable::getExN(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_n);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExN(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_n[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExN(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_n[i][j] = Ex;
  }
}


// proton
int NucleusTable::getExBinP(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_p[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinP(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_p[i] = j;
  }
}
double NucleusTable::getTotalPopP(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_p =  (const double**)((_p->second).pop_p);
    double total_pop_p=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_p += pop_p[i][j];
    }
    return total_pop_p;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopP(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_p);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopP(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_p[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopP(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_p[i][j] = pop;
  }
}
const double** NucleusTable::getExP(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_p);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExP(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_p[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExP(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_p[i][j] = Ex;
  }
}

// alpha
int NucleusTable::getExBinA(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_a[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinA(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_a[i] = j;
  }
}
double NucleusTable::getTotalPopA(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_a =  (const double**)((_p->second).pop_a);
    double total_pop_a=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_a += pop_a[i][j];
    }
    return total_pop_a;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopA(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_a);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopA(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_a[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopA(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_a[i][j] = pop;
  }
}
const double** NucleusTable::getExA(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_a);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExA(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_a[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExA(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_a[i][j] = Ex;
  }
}

// deuteron
int NucleusTable::getExBinD(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_d[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinD(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_d[i] = j;
  }
}
double NucleusTable::getTotalPopD(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_d =  (const double**)((_p->second).pop_d);
    double total_pop_d=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_d += pop_d[i][j];
    }
    return total_pop_d;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopD(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_d);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopD(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_d[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopD(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_d[i][j] = pop;
  }
}
const double** NucleusTable::getExD(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_d);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExD(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_d[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExD(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_d[i][j] = Ex;
  }
}

// triton
int NucleusTable::getExBinT(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_t[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinT(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_t[i] = j;
  }
}
double NucleusTable::getTotalPopT(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_t =  (const double**)((_p->second).pop_t);
    double total_pop_t=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_t += pop_t[i][j];
    }
    return total_pop_t;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopT(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_t);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopT(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_t[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopT(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_t[i][j] = pop;
  }
}
const double** NucleusTable::getExT(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_t);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExT(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_t[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExT(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_t[i][j] = Ex;
  }
}


// he3
int NucleusTable::getExBinH(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (int) ((_p->second).Ex_bin_h[i]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExBinH(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_bin_h[i] = j;
  }
}

double NucleusTable::getTotalPopH(const char* name, int i){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    const double** pop_h =  (const double**)((_p->second).pop_h);
    double total_pop_h=0.;
    for(int j=0;j<(int)(_p->second).array;j++){
      total_pop_h += pop_h[i][j];
    }
    return total_pop_h;
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
const double** NucleusTable::getPopH(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).pop_h);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getPopH(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).pop_h[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setPopH(const char* name, int i, int j, double pop){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).pop_h[i][j] = pop;
  }
}
const double** NucleusTable::getExH(const char* name){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (const double**) ((_p->second).Ex_h);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
double NucleusTable::getExH(const char* name, int i, int j){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    return (double) ((_p->second).Ex_h[i][j]);
  }else{
    std::cerr << "ERROR : Cannnot find such nucleus " << name << std::endl;
    abort();
  }
}
void NucleusTable::setExH(const char* name, int i, int j, double Ex){
  _p = _nucleus_table.find(name);
  if(_p!=_nucleus_table.end()){
    (_p->second).Ex_h[i][j] = Ex;
  }
}
