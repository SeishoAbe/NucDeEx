#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>

#include "NucDeExConsts.hh"
#include "NucDeExUtils.hh"
#include "NucDeExNucleusTable.hh"


///////////////
NucDeExNucleusTable::NucDeExNucleusTable()
///////////////
{
  flag_read=0;
  num_of_nuc=-1;
  //NucDeEx::Utils::SetPATH();
}

///////////////
bool NucDeExNucleusTable::ReadTables(const bool init_flag)
///////////////
{
  if(flag_read) return 0; // already read -> END

  // --- Read nucleus / separation energy tables ---//
  std::string filename1= NucDeEx::Utils::NUCDEEX_ROOT+(std::string)"/tables/nucleus/nucleus.txt";
  std::ifstream ifs(filename1);
  if(!ifs.is_open()){
    std::cerr << "ERROR : Cannot open " << filename1 << std::endl;
    //return 0;
    exit(1);
  }else{
    if(NucDeEx::Utils::fVerbose>0) std::cout << "Read: " << filename1 << std::endl;
  }
// at first get num of nucleus in the table
  int index=0;
  char buf[500];
  while(ifs.getline(buf,sizeof(buf))){
    if(buf[0]=='#') continue;
    index++;
  }
  num_of_nuc = index;

  // create array for nucleus
  _nucleus = new NucDeExNucleus[num_of_nuc];
  for(int i=0;i<num_of_nuc;i++){
    _nucleus[i].Init(init_flag);
  }

  // then read again
  char Name[5];
  int z,n;
  int maxlevelsbin;
  index=0;
  ifs.clear();
  ifs.seekg(0);
  while(ifs.getline(buf,sizeof(buf))){
    if(buf[0]=='#') continue;
    std::istringstream(buf) >> Name >> z >> n >> maxlevelsbin;
    _nucleus_id.insert(std::make_pair(Name,index));
    //_nucleus[index].name = Name;
    strcpy(_nucleus[index].name,Name);
    _nucleus[index].Z = z;
    _nucleus[index].N = n;
    _nucleus[index].A = z+n;
    _nucleus[index].id = index;
    _nucleus[index].maxlevelsbin = maxlevelsbin;

    // read separation energy table
    std::string filename2= NucDeEx::Utils::NUCDEEX_ROOT
                      + (std::string)"/tables/separation_energy/separation_energy_"
                      + (std::string)Name
                      + (std::string)".txt";
    std::ifstream ifs2(filename2);
    if(!ifs2.is_open()){
      std::cerr << "Warning: We do not have separation energy file for " << Name << std::endl;
    }else{
      if(NucDeEx::Utils::fVerbose>0) std::cout << "Read: " << filename2 << std::endl;
      _nucleus[index].flag_s = 1;
      int particle_index=0;
      while(ifs2.getline(buf,sizeof(buf))){
        if(buf[0]=='#') continue;
        std::istringstream(buf) >> _nucleus[index].S[particle_index];
        if(NucDeEx::Utils::fVerbose>0) std::cout << "SE: " << NucDeEx::particle_name[particle_index].substr(0,1) << " = " << _nucleus[index].S[particle_index] << std::endl;
        particle_index++;
      }
      ifs2.close();
    }
    index++;
  }
  ifs.close();
  flag_read=1;

  return true;
}

///////////////
NucDeExNucleus* NucDeExNucleusTable::GetNucleusPtr(int id)
///////////////
{
  if(id<0) return NULL;
  NucDeExNucleus* ptr = _nucleus;
  return ptr+id;
}

///////////////
NucDeExNucleus* NucDeExNucleusTable::GetNucleusPtr(const char* name)
///////////////
{
  return GetNucleusPtr(getID(name));
}

///////////////
NucDeExNucleus* NucDeExNucleusTable::GetNucleusPtr(int Z, int N)
///////////////
{
  if((unsigned)Z>=sizeof(nuc_name)/sizeof(char*) || Z<0) return NULL;
  //std::cout << "SIZE " << sizeof(nuc_name)/sizeof(char*) << std::endl;
  std::ostringstream os;
  os << Z+N << nuc_name[Z];
  return GetNucleusPtr(os.str().c_str());
}

///////////////
NucDeExNucleus* NucDeExNucleusTable::GetNucleusPtrPDG(int PDG)
///////////////
{
  int A = (PDG%1000)/10;
  int Z = ((PDG%1000000)-A*10)/10000;
  return GetNucleusPtr(Z,A-Z);
}

///////////////
int NucDeExNucleusTable::getID(const char* name)
///////////////
{
  _p_id = _nucleus_id.find(name);
  if(_p_id!=_nucleus_id.end()){
    return (int) ((_p_id->second));
  }else{
    if(NucDeEx::Utils::fVerbose>0) std::cerr << "Warning @ NucDeExNucleusTable: Cannot find such nucleus " << name << std::endl;
    return -1;
  }
}
