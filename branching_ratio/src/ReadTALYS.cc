#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <algorithm>

#include "Nucleus.hh"
#include "NucleusTable.hh"
#include "ReadTALYS.hh"

using namespace std;

///////////////////////////
ReadTALYS::ReadTALYS(const char* filename, NucleusTable* nuc)
///////////////////////////
{
	_filename = filename;
	_nucleus_table = nuc;
}

///////////////////////////
void ReadTALYS::SetKeywords()
///////////////////////////
{
	keyword_multiple_emission = "MULTIPLE EMISSION";
	keyword_population = "Population of Z=";
	keyword_N="N=";
	keyword_before_decay="before decay:";
}

///////////////////////////
bool ReadTALYS::Read()
///////////////////////////
{
	SetKeywords();
	_ifs = new ifstream(_filename);
  if(!_ifs->is_open()){
		cerr << "ERROR: Cannot open " << _filename << endl;
		return 0;
	}
	bool status=0;
	if(ReadPopulation() && ReadDecay()) status = 1;
	else status=1;
	_ifs->close();
	delete _ifs;
	return status;
}

///////////////////////////
bool ReadTALYS::ReadPopulation()
///////////////////////////
{
	cout << "ReadTALYS::ReadPopulation()" << endl;
  char buf[500];
	// skip untill ...
  while(_ifs->getline(buf,sizeof(buf))){
		if(strstr(buf,keyword_multiple_emission.c_str())!=NULL) break;
    else continue; 
	}
	//cout << buf << endl;

	// read continued sentences
  while(_ifs->getline(buf,sizeof(buf))){
		char* tmp=strstr(buf,keyword_population.c_str());

		// --- Population info (starting from 'keyword_population') is found
		// ->  get total population
		if(tmp!=NULL){
			string st = string(tmp);
			st = st.substr(keyword_population.length()); // remove the keyword
			int z, n;
			double total_pop;
			//cout << endl;
			//cout << tmp << endl;
			//cout << st << endl;
			istringstream(st) >> z; // obtain z
			//cout << st << endl;
			int find_N = st.find(keyword_N.c_str());
			if(find_N == string::npos){
				cerr << "something wrong in finding " << keyword_N << endl;
				return 0;
			}
			st = st.substr(find_N+keyword_N.length());
			//cout << st << endl;
			istringstream(st) >> n; // obtain n

			int find_name_s = st.find("(");
			int find_name_e = st.find(")");
			if(find_name_s == string::npos||find_name_e==string::npos){
				cerr << "something wrong in finding name" << endl;
				return 0;
			}
			string name = st.substr(find_name_s+1,find_name_e-find_name_s-1); // obtain name
			name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end()); // remove space

			int find_before_decay = st.find(keyword_before_decay.c_str());
			if(find_before_decay == string::npos){
				cerr << "something wrong in finding " << keyword_before_decay << endl;
				return 0;
			}
			st = st.substr(find_before_decay+keyword_before_decay.length());
			istringstream(st) >> total_pop;

			// Fill params
			//int ID = _nucleus_table->getID(name.c_str());
			Nucleus* nuc = _nucleus_table->GetNucleusPtr(name.c_str());
			nuc->total_pop = total_pop;
			cout <<  nuc->Z << " " << nuc->N << " " << nuc->A << " " << nuc->total_pop << endl;
		}



	}
	
	return 1;
}
///////////////////////////
bool ReadTALYS::ReadDecay()
///////////////////////////
{
	return 1;
}
