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
#include "consts.hh"

using namespace std;

///////////////////////////
ReadTALYS::ReadTALYS(const char* filename, NucleusTable* nuc)
///////////////////////////
{
	_filename = filename;
	_nucleus_table = nuc;
	os = new ostringstream;
	_verbose=0;
}

///////////////////////////
void ReadTALYS::SetKeywords()
///////////////////////////
{
	keyword_population = new string("Population of Z=");
	keyword_N = new string("N=");
	keyword_parity = new string("Parity=");
	keyword_before_decay = new string("before decay:");
	keyword_decay = new string("Decay of Z=");
	keyword_total = new string("Total:");
	keyword_bin_mother = new string("Bin=");
	keyword_parity_mother = new string("P=");
	keyword_parity_daughter = new string("P=");
	keyword_discrete = new string("Discrete levels of Z=");
	keyword_discrete_br = new string("--->");
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

	cout << "ReadTALYS::Read()" << endl;
  char buf[500];

	int flag_mode=-1;
		// 0 -> population mode
		// 1 -> decay mode
		// 2 -> discrete level mode

	Nucleus* nuc;
	int parity_array=0, parity_array_daughter=0;
	float pop_r[bins], Ex_r[bins];
	int max_bin_r;
	float pop_total_decay; // pop for specific mother -> daughter
	bool flag_first_population=1;
  while(_ifs->getline(buf,sizeof(buf))){
		string st = string(buf);

		// --- Find discrete level info 
		int find_discrete = st.find(keyword_discrete->c_str());
		if(find_discrete != string::npos){
			flag_mode=2;
			if(_verbose>0){
				cout << "### Discrete mode ###" << endl;
			}
			st = st.substr(find_discrete+keyword_discrete->length());
			int z,n;
			istringstream(st) >> z;
			int find_N = st.find(keyword_N->c_str());
			if(find_N == string::npos){
				cerr << "something wrong in finding " << keyword_N << endl;
				return 0;
			}
			st = st.substr(find_N+keyword_N->length());
			istringstream(st) >> n; // obtain n
			nuc = _nucleus_table->GetNucleusPtr(z, n);
			if(_verbose>0){
				cout << nuc->name << " " << z << " " << n << endl;
			}
		}else if(flag_mode==2){// discrete mode
			int find_discrete_br = st.find(keyword_discrete_br->c_str());

			int d_bin;
			float d_energy, d_spin;
			string d_parity;
			if(find_discrete_br == string::npos){ // level info
				if(st.find("Number") != string::npos) continue;
				if(istringstream(st) >> d_bin >> d_energy >> d_spin >> d_parity){
					if(nuc->level_Ex_bin<d_bin) nuc->level_Ex_bin = d_bin;
					nuc->level_Ex[d_bin]=d_energy;
				}
			}else{ // br info
				st = st.substr(find_discrete_br+keyword_discrete_br->length());
				int d_bin_daughter;
				float tmp_br; // br (%)
				istringstream(st) >> d_bin_daughter >> tmp_br;
				nuc->level_br[d_bin][d_bin_daughter] = tmp_br;
				if(_verbose>0){
					cout << d_bin << " ---> " << d_bin_daughter  << " : " << nuc->level_br[d_bin][d_bin_daughter] << endl;
				}
			}
		}
		
		// --- Find total population info (starting from 'keyword_population')
		// ->  get total population
		int find_population = st.find(keyword_population->c_str());
		if(find_population != string::npos){
			if(_verbose>0 && flag_mode!=0){
				cout << "### Population mode ###" << endl;
			}
			flag_mode=0;
			st = st.substr(find_population+keyword_population->length()); // remove the keyword
			int z, n;
			float total_pop;
			istringstream(st) >> z; // obtain z
			int find_N = st.find(keyword_N->c_str());
			if(find_N == string::npos){
				cerr << "something wrong in finding " << keyword_N << endl;
				return 0;
			}
			st = st.substr(find_N+keyword_N->length());
			istringstream(st) >> n; // obtain n

			nuc = _nucleus_table->GetNucleusPtr(z, n);
			nuc->flag_data = 1; // this nucleus has population data
			if(flag_first_population==1){ // first population -> target nucleus
				nuc->flag_target=1;
				flag_first_population=0;//turn off
			}
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
				nuc->sum_pop = total_pop;
				if(_verbose>0){
					cout <<  "total population: " << nuc->name << " " << nuc->Z << " " << nuc->N << " " 
							 << nuc->A << " " << nuc->sum_pop << endl;
				}
			}else{ // parity dependent
				st = st.substr(find_parity+keyword_parity->length());
				int parity;
				istringstream(st) >> parity;
				if(parity<0) parity_array=0;
				else parity_array=1;

				int find_before_decay = st.find(keyword_before_decay->c_str());
				if(find_before_decay == string::npos){
					cerr << "something wrong in finding " << keyword_before_decay << endl;
					return 0;
				}
				float pop;
				st = st.substr(find_before_decay+keyword_before_decay->length());
				istringstream(st) >> pop;
				nuc->total_pop[parity_array] = pop;
				if(_verbose>0){
					cout <<  "parity dependent population: " << nuc->name << " (" << nuc->Z << "," << nuc->N 
							 << "," << nuc->A << ") " << nuc->total_pop[parity_array] << endl;
				}
			}
			//line_population=0;
		}else{ // cannot find total population 
			int bin_mother, parity_mother, parity_daughter, daughter_id;
			int find_decay = st.find(keyword_decay->c_str());
			int find_total = st.find(keyword_total->c_str());
			if(find_decay == string::npos && flag_mode==0){
				// after total population info (polulation mode), but it is not decay info.
				// -> excitation energy & pop info.
				int bin;
				float Ex,pop, pop_p,pop_spin[9]; 
				// excitation enregy and pop info is found
				// pop : parity sum popluation. ( =  pop_parity_negative + pop_parity_positive )
				// pop_p : parity dependenet ( = sum(pop_spin[] )
				if(istringstream(st) >> bin >> Ex >> pop >> pop_p >> pop_spin[0] >> pop_spin[1] >> pop_spin[2]
							>> pop_spin[3] >> pop_spin[4] >> pop_spin[5] >> pop_spin[6]
							>> pop_spin[7] >> pop_spin[8]){
					// this is the first time we found it..
					// Usually it is negative parity.
					// the distribution is the same as positive parity
					if(nuc->Ex[parity_array][bin]<0.){
						nuc->Ex_bin[parity_array] = bin+1;
						nuc->Ex[parity_array][bin]=Ex;
						nuc->pop[parity_array][bin]=pop_p;
						if(_verbose>1){
							cout << "Parity(" << parity_array << ")" << ": excitation & pop: " << nuc->name << " " << nuc->Ex_bin[parity_array] 
									 << " " << setw(9) << nuc->Ex[parity_array][bin] << " " << setw(9) << nuc->pop[parity_array][bin]  << endl;
						}
					}
				}
			}else if (find_decay!=string::npos){ // decay info was found
				//if(_verbose>0){
				//	cout << "### Decay mode ###" << endl;
				//}
				flag_mode=1;
				int find_bin_mother = st.find(keyword_bin_mother->c_str());
				if(find_bin_mother == string::npos){
					cerr << "unexpected behavior" << endl;
					return 0;
				}
				st = st.substr(find_bin_mother+keyword_bin_mother->length());
				istringstream(st) >> bin_mother; // obtain mother's bin

				nuc->flag_decay_data[bin_mother]=1; // this nucleus have decay data!

				int find_parity_mother = st.find(keyword_parity_mother->c_str());
				st = st.substr(find_parity_mother+keyword_parity_mother->length());
				istringstream(st) >> parity_mother;
				if(parity_mother<0) parity_array=0;
				else parity_array=1;

				int find_parity_daughter = st.find(keyword_parity_daughter->c_str());
				st = st.substr(find_parity_daughter+keyword_parity_daughter->length());
				istringstream(st) >> parity_daughter;
				if(parity_daughter<0) parity_array_daughter=0;
				else parity_array_daughter=1;

				for(int i=0;i<num_particle;i++){
					if(st.find(particle_name[i].c_str())!=string::npos){
						daughter_id=i; // obtain daugher info
						break;
					}
				}
			}else if(flag_mode==1 && find_total!=string::npos){ // find total info just after decay mode
				st = st.substr(find_total+keyword_total->length());
				istringstream(st) >> pop_total_decay; // obtain pop for the decay
						// this parameter wil be sum of transition pop for both parity (negative + positive)
				if(_verbose>1){
					cout << "Total pop (decay): " << pop_total_decay << endl;
				}
			}else if(flag_mode==1){ // not decay info, but decay mode -> decay pop info
				int bin;
				float Ex,pop_spin[10]; 
				if(istringstream(st) >> bin >> Ex >> pop_spin[0] >> pop_spin[1] >> pop_spin[2]
							>> pop_spin[3] >> pop_spin[4] >> pop_spin[5] >> pop_spin[6]
							>> pop_spin[7] >> pop_spin[8] >> pop_spin[9]){

					// reject junk sentense
					int junk;
					float junkf;
					if(istringstream(st) >> junk >> junkf >> junkf >> junk >> junkf >> junkf
							>> junkf >> junkf >> junkf >> junkf >> junkf >> junkf >> junkf >> junkf) continue;
		
					// init if this is the first time to read...
					for(int i=0;i<bins && bin==0 && parity_array_daughter==0;i++){
						pop_r[i] = 0;
						Ex_r[i] =0;
						max_bin_r = 0;
					}

					// store values
					Ex_r[bin] =Ex;
					for(int i=0;i<10;i++){
						pop_r[bin] += pop_spin[i];
					}
					// store max bin if negative parity
					if(parity_array_daughter==0 && max_bin_r<bin) max_bin_r = bin;
					// fill values if the end of positve parity (last!)
					if(parity_array_daughter==1 && max_bin_r==bin){
						nuc->Ex_bin_p[daughter_id][bin_mother] = max_bin_r+1;
						float pop_total_decay_r=0;
						//for(int i=0;i<=max_bin_r;i++){
						for(int i=0;i<bins;i++){
							nuc->pop_p[daughter_id][bin_mother][i] = pop_r[i];
							nuc->Ex_p[daughter_id][bin_mother][i] = Ex_r[i];
							pop_total_decay_r += pop_r[i];
						}

						// check total pop for 
						if( (pop_total_decay_r==0 && pop_total_decay!=0)
								|| (pop_total_decay_r!=0 && abs(pop_total_decay-pop_total_decay_r)/pop_total_decay_r>check_criteria)){
							cerr << "WARNING: total pop for decay (transition pop) is not reproduced!" << endl;
							cerr << "Total pop (decay): " << pop_total_decay << endl;
							cerr << "Check total pop (decay): " << pop_total_decay_r << endl;
							cerr << nuc->name << " " << bin_mother << " " << particle_name[daughter_id] << endl;
							//return 0; // this sometimes happens.. due to bugs in TALYS
						}
						flag_mode=-1; // turn off decay mode
					}
				}
			}
		}// end of "cannot find total population info"
	} // end of getline loop
	_ifs->close();
	delete _ifs;
	return 1;
}
