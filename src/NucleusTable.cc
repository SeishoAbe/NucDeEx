#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>

#include "NucleusTable.hh"
#include "consts.hh"

#ifdef INCL_DEEXCITATION_NUCDEEX
#include "G4INCLConfig.hh"
#endif

///////////////
NucleusTable::NucleusTable()
///////////////
{
  verbose=1;
  num_of_nuc=-1;
  const char* env = getenv("NUCDEEX_TABLES");
  if(env!=NULL){
    PATH_NucDeEx_table = env;
  }else{
    std::cerr << "PATH to nucleus table is not specified" << std::endl;
    exit(1);
  }
}

#ifdef INCL_DEEXCITATION_NUCDEEX
///////////////
NucleusTable::NucleusTable(G4INCL::Config *config)
///////////////
{
  verbose=1;
  num_of_nuc=-1;
  PATH_NucDeEx_table = config->getNucDeExDataFilePath() + "/tables";
  //std::cout << PATH_NucDeEx_table.c_str() << std::endl;
}
#endif

///////////////
bool NucleusTable::ReadTables(const bool init_flag)
///////////////
{
	// --- Read nucleus / separation energy tables ---//
	std::string filename1= PATH_NucDeEx_table+(std::string)"/nucleus/nucleus.txt";
  std::ifstream ifs(filename1);
  if(!ifs.is_open()){
    std::cerr << "ERROR : Cannot open " << filename1 << std::endl;
		return 0;
  }else{
		if(verbose>0) std::cout << "Read: " << filename1 << std::endl;
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
	_nucleus = new Nucleus[num_of_nuc];
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
		std::string filename2= PATH_NucDeEx_table
                      + (std::string)"/separation_energy/separation_energy_"
                      + (std::string)Name
                      + (std::string)".txt";
		std::ifstream ifs2(filename2);
		if(!ifs2.is_open()){
			std::cerr << "Warning: We do not have separation energy file for " << Name << std::endl;
		}else{
			if(verbose>0) std::cout << "Read: " << filename2 << std::endl;
			_nucleus[index].flag_s = 1;
			int particle_index=0;
			while(ifs2.getline(buf,sizeof(buf))){
				if(buf[0]=='#') continue;
				std::istringstream(buf) >> _nucleus[index].S[particle_index];
				if(verbose>0) std::cout << "SE: " << NucDeEx::particle_name[particle_index].substr(0,1) << " = " << _nucleus[index].S[particle_index] << std::endl;
				particle_index++;
			}
			ifs2.close();
		}
		index++;
  }
	ifs.close();

  return true;
}

///////////////
Nucleus* NucleusTable::GetNucleusPtr(int id)
///////////////
{
	if(id<0) return NULL;
	Nucleus* ptr = _nucleus;
	return ptr+id;
}

///////////////
Nucleus* NucleusTable::GetNucleusPtr(const char* name)
///////////////
{
	return GetNucleusPtr(getID(name));
}

///////////////
Nucleus* NucleusTable::GetNucleusPtr(int Z, int N)
///////////////
{
	if((unsigned)Z>=sizeof(nuc_name)/sizeof(char*) || Z<0) return NULL;
	//std::cout << "SIZE " << sizeof(nuc_name)/sizeof(char*) << std::endl;
	std::ostringstream os;
	os << Z+N << nuc_name[Z];
	return GetNucleusPtr(os.str().c_str());
}

///////////////
Nucleus* NucleusTable::GetNucleusPtrPDG(int PDG)
///////////////
{
	int A = (PDG%1000)/10;
	int Z = ((PDG%1000000)-A*10)/10000;
	return GetNucleusPtr(Z,A-Z);
}

///////////////
int NucleusTable::getID(const char* name)
///////////////
{
  _p_id = _nucleus_id.find(name);
  if(_p_id!=_nucleus_id.end()){
    return (int) ((_p_id->second));
  }else{
    if(verbose>0) std::cerr << "Warning @ NucleusTable: Cannot find such nucleus " << name << std::endl;
		return -1;
  }
}
