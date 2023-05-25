#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>

#include "NucleusTable.hh"
using namespace std;

///////////////
bool NucleusTable::ReadTables()
///////////////
{
	// --- Read nucleus / separation energy tables ---//
	string filename1= (string)getenv("TALYS_WORK_TABLES")+(string)"/nucleus/nucleus.txt";
  ifstream ifs(filename1);
  if(!ifs.is_open()){
    cerr << "ERROR : Cannot open " << filename1 << endl;
		return 0;
  }else{
		cout << "Read: " << filename1 << endl;
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

	// then read again
  char Name[5];
  int z,n;
	index=0;
  ifs.clear();
  ifs.seekg(0);
  while(ifs.getline(buf,sizeof(buf))){
    if(buf[0]=='#') continue;
    istringstream(buf) >> Name >> z >> n;
		_nucleus_id.insert(make_pair(Name,index));
		//_nucleus[index].name = Name;
		strcpy(_nucleus[index].name,Name);
		_nucleus[index].Z = z;
		_nucleus[index].N = n;
		_nucleus[index].A = z+n;
		_nucleus[index].id = index;

		// read separation energy table
		string filename2= (string)getenv("TALYS_WORK_TABLES")
														+ (string)"/separation_energy/separation_energy_"
														+ (string)Name
														+ (string)".txt";
		ifstream ifs2(filename2);
		if(!ifs2.is_open()){
			cerr << "We do not have separation energy file: " << filename2 << endl;
		}else{
			cout << "Read: " << filename2 << endl;
			_nucleus[index].flag_s = 1;
			while(ifs2.getline(buf,sizeof(buf))){
				if(buf[0]=='#') continue;
				istringstream(buf) >> _nucleus[index].S[0] >> _nucleus[index].S[1]
					>> _nucleus[index].S[2] >> _nucleus[index].S[4] >> _nucleus[index].S[5]
					>> _nucleus[index].S[6] >> _nucleus[index].S[3];
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
int NucleusTable::getID(const char* name)
///////////////
{
  _p_id = _nucleus_id.find(name);
  if(_p_id!=_nucleus_id.end()){
    return (int) ((_p_id->second));
  }else{
    cerr << "ERROR : Cannnot find such nucleus " << name << endl;
    abort();
  }
}
