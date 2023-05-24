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
	os = new ostringstream;
}

///////////////////////////
void ReadTALYS::SetKeywords()
///////////////////////////
{
	keyword_multiple_emission = new string("MULTIPLE EMISSION");
	keyword_population = new string("Population of Z=");
	keyword_N = new string("N=");
	keyword_parity = new string("Parity=");
	keyword_before_decay = new string("before decay:");
	keyword_decay = new string("Decay of Z=");
	keyword_bin_mother = new string("Bin=");
	
	particle_name = new string[7];
	particle_name[0]= "gamma";
	particle_name[1]= "neutron";
	particle_name[2]= "proton";
	particle_name[3]= "alpha";
	particle_name[4]= "deuteron";
	particle_name[5]= "triton";
	particle_name[6]= "helium-3";
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
		if(strstr(buf,keyword_multiple_emission->c_str())!=NULL) break;
    else continue; 
	}

	// read continued sentences
	int line_population=0;
	Nucleus* nuc;
  while(_ifs->getline(buf,sizeof(buf))){
		string st = string(buf);
		
		// --- Find total population info (starting from 'keyword_population')
		// ->  get total population
		int find_population = st.find(keyword_population->c_str());
		if(find_population != string::npos){
			st = st.substr(find_population+keyword_population->length()); // remove the keyword
			int z, n;
			float total_pop;
			//float pop[2]={0};
			istringstream(st) >> z; // obtain z
			int find_N = st.find(keyword_N->c_str());
			if(find_N == string::npos){
				cerr << "something wrong in finding " << keyword_N << endl;
				return 0;
			}
			st = st.substr(find_N+keyword_N->length());
			istringstream(st) >> n; // obtain n

			int find_name_s = st.find("(");
			int find_name_e = st.find(")");
			if(find_name_s == string::npos||find_name_e==string::npos){
				cerr << "something wrong in finding name" << endl;
				return 0;
			}
			string name = st.substr(find_name_s+1,find_name_e-find_name_s-1); // obtain name
			name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end()); // remove space

			nuc = _nucleus_table->GetNucleusPtr(name.c_str());
			int find_parity= st.find(keyword_parity->c_str());
			if(find_parity == string::npos){ // parity plus & minus
				int find_before_decay = st.find(keyword_before_decay->c_str());
				if(find_before_decay == string::npos){
					cerr << "something wrong in finding " << keyword_before_decay << endl;
					return 0;
				}
				st = st.substr(find_before_decay+keyword_before_decay->length());
				istringstream(st) >> total_pop;

				// Fill params
				nuc->total_pop = total_pop;
				cout <<  "total population: " << nuc->name << " " << nuc->Z << " " << nuc->N << " " << nuc->A << " " << nuc->total_pop << endl;
			}else{ // parity dependent
			/*
				st = st.substr(find_parity+keyword_parity->length());
				int parity;
				istringstream(st) >> parity;
				int array=0;
				if(parity<0) array=0;
				else array=1;

				int find_before_decay = st.find(keyword_before_decay->c_str());
				if(find_before_decay == string::npos){
					cerr << "something wrong in finding " << keyword_before_decay << endl;
					return 0;
				}
				st = st.substr(find_before_decay+keyword_before_decay->length());
				istringstream(st) >> pop[array];
				nuc->pop[array] = pop[array];
				cout <<  nuc->name << " " << nuc->Z << " " << nuc->N << " " << nuc->A << " " << nuc->pop[array] << endl;
			*/
			}
			line_population=0;
		}else{ // cannot find total population info
			line_population++; // count line after total pop info
			if(line_population<skip_after_population) continue;

			int find_decay = st.find(keyword_decay->c_str());
			if(find_decay == string::npos){
				// after total population info, but it is not decay info.
				// -> excitation energy & pop info.
				string st = buf;
				int bin;
				float Ex,pop, pop_p,pop_spin[9]; // pop_p = sum(pop_spin[i]), pop = sum(pop_p), parity sum
				// excitation enregy and pop info is found
				if(istringstream(st) >> bin >> Ex >> pop >> pop_spin[0] >> pop_spin[1] >> pop_spin[2]
							>> pop_spin[3] >> pop_spin[4] >> pop_spin[5] >> pop_spin[6]
							>> pop_spin[7] >> pop_spin[8] >> pop_spin[9]){
					// this is the first time we found it..
					// Usually it is negative parity.
					// the distribution is the same as positive parity
					if(nuc->Ex[bin]<0.){
						nuc->Ex_bin = bin;
						nuc->Ex[bin]=Ex;
						nuc->pop[bin]=pop;
						cout << "excitationa and pop: " << nuc->name << " " << nuc->Ex_bin << " " << setw(9) << nuc->Ex[bin]
								 << " " << setw(9) << nuc->pop[bin] << endl;
					}
				}
			}else{ // It is decay info
				int find_bin_mother = st.find(keyword_bin_mother->c_str());
				if(find_bin_mother == string::npos){
					cerr << "unexpected behavior" << endl;
					return 0;
				}
				st = st.substr(find_bin_mother+keyword_bin_mother->length());
				int bin_mother;
				istringstream(st) >> bin_mother;

				int daughter_id;
				for(int i=0;i<num_particle;i++){
					if(st.find(particle_name[i].c_str())!=string::npos){
						daughter_id=i;
						break;
					}
				}

			}

		}


		// --- Find 

	} 
	return 1;
}
///////////////////////////
bool ReadTALYS::ReadDecay()
///////////////////////////
{
	return 1;
}
