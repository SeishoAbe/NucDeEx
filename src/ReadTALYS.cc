#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <algorithm>

#include "NucDeExNucleus.hh"
#include "ReadTALYS.hh"

///////////////////////////
ReadTALYS::ReadTALYS(const char* filename, NucDeExNucleusTable* nuc)
///////////////////////////
{
  _filename = filename;
  _nucleus_table = nuc;
  os = new std::ostringstream;
  _verbose=0;
}

///////////////////////////
void ReadTALYS::SetKeywords()
///////////////////////////
{
  keyword_population = new std::string("Population of Z=");
  keyword_N = new std::string("N=");
  keyword_parity = new std::string("Parity=");
  keyword_before_decay = new std::string("before decay:");
  keyword_decay = new std::string("Decay of Z=");
  keyword_total = new std::string("Total:");
  keyword_bin_mother = new std::string("Bin=");
  keyword_parity_mother = new std::string("P=");
  keyword_parity_daughter = new std::string("P=");
  keyword_discrete = new std::string("Discrete levels of Z=");
  keyword_discrete_br = new std::string("--->");
}

///////////////////////////
bool ReadTALYS::Read()
///////////////////////////
{
  SetKeywords();
  _ifs = new std::ifstream(_filename);
  if(!_ifs->is_open()){
    std::cerr << "ERROR: Cannot open " << _filename << std::endl;
    return 0;
  }

  std::cout << "ReadTALYS::Read()" << std::endl;
  char buf[500];

  int flag_mode=-1;
    // 0 -> population mode
    // 1 -> decay mode
    // 2 -> discrete level mode

  NucDeExNucleus* nuc;
  int parity_array=0, parity_array_daughter=0;
  float pop_r[NucDeEx::bins], Ex_r[NucDeEx::bins];
  int max_bin_r;
  float pop_total_decay; // pop for specific mother -> daughter
  bool flag_first_population=1;
  while(_ifs->getline(buf,sizeof(buf))){
    std::string st = std::string(buf);

    // --- Find discrete level info 
    int find_discrete = st.find(keyword_discrete->c_str());
    if(find_discrete != std::string::npos){
      flag_mode=2;
      if(_verbose>0){
        std::cout << "### Discrete mode ###" << std::endl;
      }
      st = st.substr(find_discrete+keyword_discrete->length());
      int z,n;
      std::istringstream(st) >> z;
      int find_N = st.find(keyword_N->c_str());
      if(find_N == std::string::npos){
        std::cerr << "something wrong in finding " << keyword_N << std::endl;
        return 0;
      }
      st = st.substr(find_N+keyword_N->length());
      std::istringstream(st) >> n; // obtain n
      nuc = _nucleus_table->GetNucleusPtr(z, n);
      if(_verbose>0){
        std::cout << nuc->name << " " << z << " " << n << std::endl;
      }
    }else if(flag_mode==2){// discrete mode
      int find_discrete_br = st.find(keyword_discrete_br->c_str());

      int d_bin;
      float d_energy, d_spin;
      std::string d_parity;
      if(find_discrete_br == std::string::npos){ // level info
        if(st.find("Number") != std::string::npos) continue;
        if(std::istringstream(st) >> d_bin >> d_energy >> d_spin >> d_parity){
          if(nuc->level_Ex_bin<d_bin) nuc->level_Ex_bin = d_bin;
          nuc->level_Ex[d_bin]=d_energy;
        }
      }else{ // br info
        st = st.substr(find_discrete_br+keyword_discrete_br->length());
        int d_bin_daughter;
        float tmp_br; // br (%)
        std::istringstream(st) >> d_bin_daughter >> tmp_br;
        nuc->level_br[d_bin][d_bin_daughter] = tmp_br;
        if(_verbose>0){
          std::cout << d_bin << " ---> " << d_bin_daughter  << " : " << nuc->level_br[d_bin][d_bin_daughter] << std::endl;
        }
      }
    }
    
    // --- Find total population info (starting from 'keyword_population')
    // ->  get total population
    int find_population = st.find(keyword_population->c_str());
    if(find_population != std::string::npos){
      if(_verbose>0 && flag_mode!=0){
        std::cout << "### Population mode ###" << std::endl;
      }
      flag_mode=0;
      st = st.substr(find_population+keyword_population->length()); // remove the keyword
      int z, n;
      float total_pop;
      std::istringstream(st) >> z; // obtain z
      int find_N = st.find(keyword_N->c_str());
      if(find_N == std::string::npos){
        std::cerr << "something wrong in finding " << keyword_N << std::endl;
        return 0;
      }
      st = st.substr(find_N+keyword_N->length());
      std::istringstream(st) >> n; // obtain n

      nuc = _nucleus_table->GetNucleusPtr(z, n);
      nuc->flag_data = 1; // this nucleus has population data
      if(flag_first_population==1){ // first population -> target nucleus
        nuc->flag_target=1;
        flag_first_population=0;//turn off
      }
      int find_parity= st.find(keyword_parity->c_str());
      if(find_parity == std::string::npos){ // parity plus & minus
        int find_before_decay = st.find(keyword_before_decay->c_str());
        if(find_before_decay == std::string::npos){
          std::cerr << "something wrong in finding " << keyword_before_decay << std::endl;
          return 0;
        }
        st = st.substr(find_before_decay+keyword_before_decay->length());
        std::istringstream(st) >> total_pop;

        // Fill params
        nuc->sum_pop = total_pop;
        if(_verbose>0){
          std::cout <<  "total population: " << nuc->name << " " << nuc->Z << " " << nuc->N << " " 
               << nuc->A << " " << nuc->sum_pop << std::endl;
        }
      }else{ // parity dependent
        st = st.substr(find_parity+keyword_parity->length());
        int parity;
        std::istringstream(st) >> parity;
        if(parity<0) parity_array=0;
        else parity_array=1;

        int find_before_decay = st.find(keyword_before_decay->c_str());
        if(find_before_decay == std::string::npos){
          std::cerr << "something wrong in finding " << keyword_before_decay << std::endl;
          return 0;
        }
        float pop;
        st = st.substr(find_before_decay+keyword_before_decay->length());
        std::istringstream(st) >> pop;
        nuc->total_pop[parity_array] = pop;
        if(_verbose>0){
          std::cout <<  "parity dependent population: " << nuc->name << " (" << nuc->Z << "," << nuc->N 
               << "," << nuc->A << ") " << nuc->total_pop[parity_array] << std::endl;
        }
      }
      //line_population=0;
    }else{ // cannot find total population 
      int bin_mother, parity_mother, parity_daughter, daughter_id;
      int find_decay = st.find(keyword_decay->c_str());
      int find_total = st.find(keyword_total->c_str());
      if(find_decay == std::string::npos && flag_mode==0){
        // after total population info (polulation mode), but it is not decay info.
        // -> excitation energy & pop info.
        int bin;
        float Ex,pop, pop_p,pop_spin[9]; 
        // excitation enregy and pop info is found
        // pop : parity sum popluation. ( =  pop_parity_negative + pop_parity_positive )
        // pop_p : parity dependenet ( = sum(pop_spin[] )
        if(std::istringstream(st) >> bin >> Ex >> pop >> pop_p >> pop_spin[0] >> pop_spin[1] >> pop_spin[2]
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
              std::cout << "Parity(" << parity_array << ")" << ": excitation & pop: " << nuc->name << " " << nuc->Ex_bin[parity_array] 
                   << " " << std::setw(9) << nuc->Ex[parity_array][bin] << " " << std::setw(9) << nuc->pop[parity_array][bin]  << std::endl;
            }
          }
        }
      }else if (find_decay!=std::string::npos){ // decay info was found
        //if(_verbose>0){
        //  std::cout << "### Decay mode ###" << std::endl;
        //}
        flag_mode=1;
        int find_bin_mother = st.find(keyword_bin_mother->c_str());
        if(find_bin_mother == std::string::npos){
          std::cerr << "unexpected behavior" << std::endl;
          return 0;
        }
        st = st.substr(find_bin_mother+keyword_bin_mother->length());
        std::istringstream(st) >> bin_mother; // obtain mother's bin

        nuc->flag_decay_data[bin_mother]=1; // this nucleus have decay data!

        int find_parity_mother = st.find(keyword_parity_mother->c_str());
        st = st.substr(find_parity_mother+keyword_parity_mother->length());
        std::istringstream(st) >> parity_mother;
        if(parity_mother<0) parity_array=0;
        else parity_array=1;

        int find_parity_daughter = st.find(keyword_parity_daughter->c_str());
        st = st.substr(find_parity_daughter+keyword_parity_daughter->length());
        std::istringstream(st) >> parity_daughter;
        if(parity_daughter<0) parity_array_daughter=0;
        else parity_array_daughter=1;

        for(int i=0;i<NucDeEx::num_particle;i++){
          if(st.find(NucDeEx::particle_name[i].c_str())!=std::string::npos){
            daughter_id=i; // obtain daugher info
            break;
          }
        }
      }else if(flag_mode==1 && find_total!=std::string::npos){ // find total info just after decay mode
        st = st.substr(find_total+keyword_total->length());
        std::istringstream(st) >> pop_total_decay; // obtain pop for the decay
            // this parameter wil be sum of transition pop for both parity (negative + positive)
        if(_verbose>1){
          std::cout << "Total pop (decay): " << pop_total_decay << std::endl;
        }
      }else if(flag_mode==1){ // not decay info, but decay mode -> decay pop info
        int bin;
        float Ex,pop_spin[10]; 
        if(std::istringstream(st) >> bin >> Ex >> pop_spin[0] >> pop_spin[1] >> pop_spin[2]
              >> pop_spin[3] >> pop_spin[4] >> pop_spin[5] >> pop_spin[6]
              >> pop_spin[7] >> pop_spin[8] >> pop_spin[9]){

          // reject junk sentense
          int junk;
          float junkf;
          if(std::istringstream(st) >> junk >> junkf >> junkf >> junk >> junkf >> junkf
              >> junkf >> junkf >> junkf >> junkf >> junkf >> junkf >> junkf >> junkf) continue;
    
          // init if this is the first time to read...
          for(int i=0;i<NucDeEx::bins && bin==0 && parity_array_daughter==0;i++){
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
            for(int i=0;i<NucDeEx::bins;i++){
              nuc->pop_p[daughter_id][bin_mother][i] = pop_r[i];
              nuc->Ex_p[daughter_id][bin_mother][i] = Ex_r[i];
              pop_total_decay_r += pop_r[i];
            }

            // check total pop for 
            if( (pop_total_decay_r==0 && pop_total_decay!=0)
                || (pop_total_decay_r!=0 && abs(pop_total_decay-pop_total_decay_r)/pop_total_decay_r>check_criteria)){
              std::cerr << "WARNING: total pop for decay (transition pop) is not reproduced!" << std::endl;
              std::cerr << "Total pop (decay): " << pop_total_decay << std::endl;
              std::cerr << "Check total pop (decay): " << pop_total_decay_r << std::endl;
              std::cerr << nuc->name << " " << bin_mother << " " << NucDeEx::particle_name[daughter_id] << std::endl;
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
