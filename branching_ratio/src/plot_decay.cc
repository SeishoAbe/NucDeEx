#include "NucleusTable.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>

int main(int argc, char* argv[]){
  // --- FIXME  --- //
  std::string parent="11C";
  const int maxlevelstar=8;
  const int maxlevelsbin_g=8;
  const int maxlevelsbin_n=10;
  const int maxlevelsbin_p=10;
  const int maxlevelsbin_a=10;
  const int maxlevelsbin_d=5;
  const int maxlevelsbin_t=5;
  const int maxlevelsbin_h=5;

  const int Ex_max=60000; // keV
  const int Ex_bin_width=100; // keV
  const double max_Ex_plot=50;
  // ---------------//
  
  std::ostringstream os;
  NucleusTable* nucleus_table = NucleusTable::getPointer();
  nucleus_table->ReadTable("nucleus.txt");


// Get a table of nucleus //

// Read output //
  const char* population = " Population of Z=";
  const int n_population = 17;
  //
  const char* decay = " Decay of Z=";
  //const int n_decay = 12;
  std::ifstream ifs;
  ifs.open("output/output");
  if(!ifs.is_open()) return 1;
  bool multiple_emission=0;
  char buf[500];

  //
  int z,n,index;
  double total_pop;
  const char* name;
  std::string name_s;
  int line=-1e9, line_decay=-1e9;
  //
  int bin, P, z_r, n_r;
  double total_pop_r;
  const char* name_r;
  std::string name_r_s;
  double Ex, J, pop;
  const int array=100;
  double Ex_r[array]={0.}, pop_r[array]={0.};
  int bin_r, max_bin_r, P_r;
  int z_rr[array]={0}, n_rr[array]={0};
  int index_rr=0;
  while(ifs.getline(buf,sizeof(buf))){
    if(strncmp(buf," ########## MULTIPLE EMISSION ",30)==0)
        multiple_emission=1;
    if(multiple_emission==0) continue;
  
    // --- Read population --- //
    // Read 1st line of population 
    if(strncmp(buf,population,n_population)==0){
      std::cout << "### POPULATION MODE###" << std::endl;
      //std::cout << buf << std::endl;
      std::string st = buf;
      //std::cout << st << std::endl;
      st = st.substr(n_population,st.size());
      std::string junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8;
      //else if(std::istringstream(st) >> z >> junk1 >> n 
      if(std::istringstream(st) >> z >> junk1 >> n 
                   >> junk2 >> junk3 >> junk4 >> junk5 >> junk6 >> total_pop){
        //std::cout << st << std::endl;
        //std::cout <<"junk5=" <<  junk5 << std::endl;
        //std::cout <<"junk6=" <<  junk6 << std::endl;
        if(!(junk5=="before" && junk6=="decay:")){
          //std::cout << junk6.size() << std::endl;
          //std::cout << junk7.size() << std::endl;
          std::cerr << "ERROR in population 2: " << st << std::endl;
          return 1;
        }
      }
      else if(std::istringstream(st) >> z >> junk1 >> n 
                   >> junk2 >> junk3 >> junk4 >> junk5 >> junk6 >> junk7 >> total_pop){
        //std::cout <<"junk5=" <<  junk5 << std::endl;
        //std::cout <<"junk6=" <<  junk6 << std::endl;
        //std::cout <<"junk7=" <<  junk7 << std::endl;
        if(!(junk6=="before" && junk7=="decay:")){
          std::cerr << "ERROR in population 3 : " << st << std::endl;
          return 1;
        }
      }
      else if(std::istringstream(st) >> z >> junk1 >> n 
                   >> junk2 >> junk3 >> junk4 >> junk5 >> junk6 >> junk7 >> junk8 >> total_pop){
        //std::cout <<"junk5=" <<  junk5 << std::endl;
        //std::cout <<"junk6=" <<  junk6 << std::endl;
        //std::cout <<"junk7=" <<  junk7 << std::endl;
        //std::cout <<"junk8=" <<  junk8 << std::endl;
        if(!(junk7=="before" && junk8=="decay:")){
          std::cerr << "ERROR in population 4 : " << st << std::endl;
          return 1;
        }
      }
      else if(std::istringstream(st) >> z >> junk1 >> n 
                   >> junk2 >> junk3 >> junk4 >> junk5 >> total_pop){
        //std::cout << st << std::endl;
        //std::cout <<"junk4=" <<  junk4 << std::endl;
        //std::cout <<"junk5=" <<  junk5 << std::endl;
        if(!(junk4=="before" && junk5=="decay:")){
          //std::cout << junk6.size() << std::endl;
          //std::cout << junk7.size() << std::endl;
          std::cerr << "ERROR in population 1: " << st << std::endl;
          return 1;
        }
      }
      else{
        std::cout << "ERROR in population 5 : no data " << std::endl;
        return 1;
      }
      name  = nucleus_table->getName(z,n);
      name_s = name;
      if(nucleus_table->getTotalPop(name)==0){ 
        //std::cout << std::endl;
        nucleus_table->setTotalPop(name,total_pop);
        z_rr[index_rr] = z;
        n_rr[index_rr] = n;
        index_rr++;
      }
      index=nucleus_table->getIndex(name);
      nucleus_table->setIndex(name,index+1);
      if(index==1) line=0;
      std::cout << name << "   " << nucleus_table->getTotalPop(name) << "   " << index << std::endl;
    }// end of reading 1st line of population

    if(index==1 && line>=4){
      std::string st = buf;
      int bin;
      //double Ex,pop;
      if(std::istringstream(st) >> bin >> Ex >> pop){
        nucleus_table->setEx(name, bin, Ex);
        nucleus_table->setPop(name, bin, pop);
        nucleus_table->setExBin(name,bin);
        //std::cout << Ex << "   " << pop << std::endl;
        //std::cout << nucleus_table->getExBin(name) << "   "
        //          << nucleus_table->getEx(name, bin) << "   " 
        //          << nucleus_table->getPop(name, bin) << std::endl;
      }
    }
    // --- End of read population --- //
    


    // --- Read decay --- //
    os.str("");
    os << resetiosflags(ios_base::floatfield);
    /*
    if(z+n>=10) 
      os << " Decay of Z=" << std::setw(3) << z << " N=" << std::setw(3) << n << " ( " << std::setw(3) << name_s << " ), Bin=" ;
    else 
      os << " Decay of Z=" << std::setw(3) << z << " N=" << std::setw(3) << n << " (  " << std::setw(3) << name_s << "), Bin=" ;
      */
    if(z+n>=10) 
      os << " Decay of Z=" << std::setw(3) << z << " N=" << std::setw(3) << n << " ( " << std::setw(4) << std::left << name_s << "), Bin=" ;
    else 
      os << " Decay of Z=" << std::setw(3) << z << " N=" << std::setw(3) << n << " (  " << std::setw(3) << std::left << name_s << "), Bin=" ;
    //std::cout << std::right;
    os << std::right;
    //std::cout << buf << std::endl;
    //std::cout << os.str().c_str()  << std::endl;
    if(index>=2 && strncmp(buf,os.str().c_str(),os.str().size())==0){
      std::cout << "### Decay MODE###" << std::endl;
      std::string st = buf;
      //std::cout << st << std::endl;
      //std::cout << os.str().c_str()  << std::endl;
      st = st.substr(os.str().size(),st.size());
      //std::cout << st << std::endl;
      std::string junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8;
      if(std::istringstream(st) >> bin >> junk1 >> Ex >> junk2 >> J){
        os.str("");
        os << resetiosflags(ios_base::floatfield);
        os << std::setw(3) << bin << " Ex=" << std::setw(8) << std::fixed << std::setprecision(3) << Ex 
           << resetiosflags(ios_base::floatfield)
            << " J=" << std::setw(4) << J << " P=";
        //std::cout << os.str().c_str() << std::endl;
        st = st.substr(os.str().size(),st.size());
        //std::cout << st << std::endl;
        if(std::istringstream(st) >> P >> junk4 >> pop){
          os.str("");
          os << resetiosflags(ios_base::floatfield);
          os << std::setw(2) << P << " Pop=" << std::setw(10) << std::scientific << std::setprecision(3) 
             << std::uppercase << pop << " to bins of Z=";
          //std::cout << os.str().c_str() << std::endl;
          st = st.substr(os.str().size(),st.size());
          if(std::istringstream(st) >> z_r >> junk1 >> n_r){
            //std::cout << st << std::endl;
            name_r  = nucleus_table->getName(z_r,n_r);
            name_r_s = name_r;
            os.str("");
            if(z_r+n_r>=10) 
              os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " ( " << std::setw(4) << std::left << name_r_s << "), P=";
            else 
              os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " (  " << std::setw(3) << std::left << name_r_s << "), P=";
            os << std::right;
            /*
            if(z_r+n_r>=10) 
              os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " ( " << std::setw(3) << name_r_s << " ), P=";
            else 
              os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " (  " << std::setw(3) << name_r_s << "), P=";
            */
            st = st.substr(os.str().size(),st.size());
            if(std::istringstream(st) >> P_r){
              os.str("");
              os << resetiosflags(ios_base::floatfield);
              os << " Decay of Z=" << std::setw(3) << z ;
              if(z+n>=10)
                os << " N=" << std::setw(3) << n << " ( " << std::setw(4) << std::left << name_s << "), Bin=";
              else 
                os << " N=" << std::setw(3) << n << " (  " << std::setw(3) << std::left << name_s << "), Bin=";
              os << std::right;
              os << std::setw(3) << bin << " Ex=" << std::setw(8) << std::fixed << std::setprecision(3) << Ex 
                 << resetiosflags(ios_base::floatfield)
                 << " J=" << std::setw(4) << J << " P=" 
                 << std::setw(2) << P << " Pop=" << std::setw(10) << std::scientific << std::setprecision(3) 
                 << std::uppercase << pop << " to bins of Z=";
              if(z_r+n_r>=10) 
                os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " ( " << std::setw(4) << std::left << name_r_s << "), P="
                   << std::setw(2) << P_r;
              else 
                os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " (  " << std::setw(3) << std::left << name_r_s << "), P="
                   << std::setw(2) << std::right << P_r;
              os << std::right;
                /*
              if(z_r+n_r>=10)                    
                os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " ( " << name_r_s << " ), P=" << std::setw(2) << P_r;
              else
                os << std::setw(3) << z_r << " N=" << std::setw(3) << n_r << " (  " << name_r_s << "), P=" << std::setw(2) << P_r;
              */
              os << "   <- CHECK";
              //std::cout << os.str().c_str() << std::endl;
              //std::cout << "READ!!!" << std::endl;
              //std::cout << resetiosflags(ios_base::floatfield);// << std::endl;
              line_decay=0;
              // initialize
              if(P_r==-1){ 
                max_bin_r=0;
                for(int i=0;i<array;i++){
                  pop_r[i] = 0.0;
                  Ex_r[i] =0.0;
                }
              }
            }
          }
          //std::cout << st << std::endl;
          //std::cout << os.str().c_str() << std::endl;
        }
      }
    }
    if(index>=2 && line_decay==2){
      std::string st = buf;
      os.str("");
      os << " Total: ";
      st = st.substr(os.str().size(),st.size());
      if(std::istringstream(st) >> total_pop_r)
        std::cout << "Total pop : " << total_pop_r << std::endl;
    }
    if(index>=2 && line_decay>=6){
      int flag_decay=-1;
      if(z==z_r && n==n_r) flag_decay=0;// gamma
      else if(z==z_r && n==n_r+1) flag_decay=1; // neutron
      else if(z==z_r+1 && n==n_r) flag_decay=2; // proton
      else if(z==z_r+2 && n==n_r+2) flag_decay=3; // alpha
      else if(z==z_r+1 && n==n_r+1) flag_decay=4; // deuteron
      else if(z==z_r+1 && n==n_r+2) flag_decay=5; // triton
      else if(z==z_r+2 && n==n_r+1) flag_decay=6; // he3
      else continue;

      std::string st = buf;
      double pop_r_tmp[10];
      double Ex_r_tmp;
      if(std::istringstream(st) >> bin_r >> Ex_r_tmp >> pop_r_tmp[0] >> pop_r_tmp[1] >> pop_r_tmp[2] >> pop_r_tmp[3] >> pop_r_tmp[4]
                                 >> pop_r_tmp[5] >> pop_r_tmp[6] >> pop_r_tmp[7] >> pop_r_tmp[8] >> pop_r_tmp[9]){
        if(P_r==-1 && max_bin_r<bin_r) max_bin_r=bin_r;
        if(P_r==1 && bin_r > max_bin_r){
          line_decay=-1e9;
          continue;
        }
        if(P_r==1 && (Ex_r_tmp != Ex_r[bin_r] )){
          line_decay=-1e9;
          continue;
        }
        //std::cout << st << std::endl;
        //std::cout << bin_r << std::endl;
        Ex_r[bin_r] = Ex_r_tmp;
        double sum_pop_r=0;
        for(int i=0;i<10;i++){
          sum_pop_r += pop_r_tmp[i];
        }
        pop_r[bin_r] += sum_pop_r;
        double total_pop_r_check=0.;
        if(P_r==1 && max_bin_r==bin_r){
          for(int i=0;i<=max_bin_r;i++){
            total_pop_r_check += pop_r[i];
            std::cout << std::setw(5) <<  i << "   "  << std::setw(10) << Ex_r[i] << "   " 
                      << std::setw(10) << pop_r[i] << std::endl;
            if(flag_decay==0){ // gamma
              nucleus_table->setExBinG(name,bin, max_bin_r);
              nucleus_table->setExG(name,bin,i,Ex_r[i]);
              nucleus_table->setPopG(name,bin,i,pop_r[i]);
            }else if(flag_decay==1){ // neutron
              nucleus_table->setExBinN(name,bin, max_bin_r);
              nucleus_table->setExN(name,bin,i,Ex_r[i]);
              nucleus_table->setPopN(name,bin,i,pop_r[i]);
            }else if(flag_decay==2){ // proton
              nucleus_table->setExBinP(name,bin, max_bin_r);
              nucleus_table->setExP(name,bin,i,Ex_r[i]);
              nucleus_table->setPopP(name,bin,i,pop_r[i]);
            }else if(flag_decay==3){ // alpha
              nucleus_table->setExBinA(name,bin, max_bin_r);
              nucleus_table->setExA(name,bin,i,Ex_r[i]);
              nucleus_table->setPopA(name,bin,i,pop_r[i]);
            }else if(flag_decay==4){ // deuteron
              nucleus_table->setExBinD(name,bin, max_bin_r);
              nucleus_table->setExD(name,bin,i,Ex_r[i]);
              nucleus_table->setPopD(name,bin,i,pop_r[i]);
            }else if(flag_decay==5){ // triton
              nucleus_table->setExBinT(name,bin, max_bin_r);
              nucleus_table->setExT(name,bin,i,Ex_r[i]);
              nucleus_table->setPopT(name,bin,i,pop_r[i]);
            }else if(flag_decay==6){ // he
              nucleus_table->setExBinH(name,bin, max_bin_r);
              nucleus_table->setExH(name,bin,i,Ex_r[i]);
              nucleus_table->setPopH(name,bin,i,pop_r[i]);
            }
          }
          std::cout << total_pop_r_check << "  ";
          if (flag_decay==0) std::cout << nucleus_table->getTotalPopG(name,bin) << "   G" << std::endl;
          else if (flag_decay==1) std::cout << nucleus_table->getTotalPopN(name,bin) << "   N" << std::endl;
          else if (flag_decay==2) std::cout << nucleus_table->getTotalPopP(name,bin) << "   P" << std::endl;
          else if (flag_decay==3) std::cout << nucleus_table->getTotalPopA(name,bin) << "   A" << std::endl;
          else if (flag_decay==4) std::cout << nucleus_table->getTotalPopD(name,bin) << "   D" << std::endl;
          else if (flag_decay==5) std::cout << nucleus_table->getTotalPopT(name,bin) << "   T" << std::endl;
          else if (flag_decay==6) std::cout << nucleus_table->getTotalPopH(name,bin) << "   H" << std::endl;
          else abort();
          std::cout << std::endl;
          if(max_bin_r !=0 && total_pop_r!=0 && (total_pop_r-total_pop_r_check)/total_pop_r>0.1){
            std::cerr << "#### UNEXPECTED BEHAVIOR " << std::endl;
            abort();
          }
        }
      }
    }
    line_decay++;
    line++;
  }
  ifs.close();
  std::cout << std::endl;

  for(int i=0;i<index_rr;i++){
    std::cout << z_rr[i] << "  " << n_rr[i] << "   " << nucleus_table->getName(z_rr[i],n_rr[i]) << std::endl;
    const char* name_target =nucleus_table->getName(z_rr[i],n_rr[i]);
    std::string target = name_target;

    // Graph target
    os.str("");
    os << "txt/Br_" << target.c_str() << "_summary.root";
    TFile* rootf = new TFile(os.str().c_str(),"RECREATE");
    
    //
    double Sg=-1, Sn=-1, Sp=-1, Sd=-1, St=-1, Sh=-1, Sa=-1;
    os.str("");
    os << "/userdata/work/seisho/TALYS/work-dir_talys_1.95/separation_energy/separation_energy_" << target.c_str() <<  ".txt";
    ifs.open(os.str().c_str());
    cout << os.str().c_str() << endl;
    if(!ifs.is_open()) return 1;
    else std::cout << os.str().c_str() << std::endl;
    while(ifs.getline(buf,sizeof(buf))){
      //std::istringstream(buf) >> Sg >> Sn >> Sp >> Sd >> St >> Sh >> Sa;
      if(!( std::istringstream(buf) >> Sg >> Sn >> Sp >> Sd >> St >> Sh >> Sa)){
        cerr << "Something wrong in separation E table " << endl;
        return 1;
      }
    }
    ifs.close();
    std::cout << "Sep  E : " << Sg << " " << Sn << " " << Sp << " " << Sd << " " 
              << St << " " << Sh << " " << Sa << std::endl;
    //
    TGraph* g_target_pop = new TGraph;
    TGraph* g_target_br_g = new TGraph;
    TGraph* g_target_br_n = new TGraph;
    TGraph* g_target_br_p = new TGraph;
    TGraph* g_target_br_a = new TGraph;
    TGraph* g_target_br_d = new TGraph;
    TGraph* g_target_br_t = new TGraph;
    TGraph* g_target_br_h = new TGraph;
    int numofexbin=nucleus_table->getExBin(target.c_str());
    TGraph* g_target_br_ex_g[numofexbin];
    TGraph* g_target_br_ex_n[numofexbin];
    TGraph* g_target_br_ex_p[numofexbin];
    TGraph* g_target_br_ex_a[numofexbin];
    TGraph* g_target_br_ex_d[numofexbin];
    TGraph* g_target_br_ex_t[numofexbin];
    TGraph* g_target_br_ex_h[numofexbin];
    double total_p  = nucleus_table->getTotalPop(target.c_str());
    double max_p=0.;
    index=0;
    for(int i=0;i<=numofexbin;i++){
      bool flag_below_S=0;
      g_target_br_ex_g[i] = new TGraph;
      g_target_br_ex_n[i] = new TGraph;
      g_target_br_ex_p[i] = new TGraph;
      g_target_br_ex_a[i] = new TGraph;
      g_target_br_ex_d[i] = new TGraph;
      g_target_br_ex_t[i] = new TGraph;
      g_target_br_ex_h[i] = new TGraph;
      //
      double e = nucleus_table->getEx(target.c_str(),i);
      double p = nucleus_table->getPop(target.c_str(),i);
      g_target_pop->SetPoint(i,e,p/total_p);
      if(max_p<p/total_p) max_p=p/total_p;
      //
      if(p==0) continue;
      double total_p_g = nucleus_table->getTotalPopG(target.c_str(),i);
      if(e<Sn && e<Sp && e<Sa && e<Sd && e<St && Sh){
        flag_below_S=1;
        total_p_g=p;
      }
      double total_p_n = nucleus_table->getTotalPopN(target.c_str(),i);
      double total_p_p = nucleus_table->getTotalPopP(target.c_str(),i);
      double total_p_a = nucleus_table->getTotalPopA(target.c_str(),i);
      double total_p_d = nucleus_table->getTotalPopD(target.c_str(),i);
      double total_p_t = nucleus_table->getTotalPopT(target.c_str(),i);
      double total_p_h = nucleus_table->getTotalPopH(target.c_str(),i);
      if(total_p_g==0. && total_p_n==0. && total_p_p==0. && total_p_a==0.
         && total_p_d==0. && total_p_t==0. && total_p_h==0.){
        if(index!=0){
          double X,Y;
          g_target_br_g->GetPoint(index-1,X,Y);
          g_target_br_g->SetPoint(index,e,Y);
          g_target_br_n->GetPoint(index-1,X,Y);
          g_target_br_n->SetPoint(index,e,Y);
          g_target_br_p->GetPoint(index-1,X,Y);
          g_target_br_p->SetPoint(index,e,Y);
          g_target_br_a->GetPoint(index-1,X,Y);
          g_target_br_a->SetPoint(index,e,Y);
          g_target_br_d->GetPoint(index-1,X,Y);
          g_target_br_d->SetPoint(index,e,Y);
          g_target_br_t->GetPoint(index-1,X,Y);
          g_target_br_t->SetPoint(index,e,Y);
          g_target_br_h->GetPoint(index-1,X,Y);
          g_target_br_h->SetPoint(index,e,Y);
        }else{
          g_target_br_g->SetPoint(index,e,1);
          g_target_br_n->SetPoint(index,e,0);
          g_target_br_p->SetPoint(index,e,0);
          g_target_br_a->SetPoint(index,e,0);
          g_target_br_d->SetPoint(index,e,0);
          g_target_br_t->SetPoint(index,e,0);
          g_target_br_h->SetPoint(index,e,0);
        }
      }else{
        g_target_br_g->SetPoint(index,e,total_p_g/p);
        g_target_br_n->SetPoint(index,e,total_p_n/p);
        g_target_br_p->SetPoint(index,e,total_p_p/p);
        g_target_br_a->SetPoint(index,e,total_p_a/p);
        g_target_br_d->SetPoint(index,e,total_p_d/p);
        g_target_br_t->SetPoint(index,e,total_p_t/p);
        g_target_br_h->SetPoint(index,e,total_p_h/p);
      }

      // -- br for level -- //
      const int numofexbin_g = nucleus_table->getExBinG(target.c_str(),i);
      const int numofexbin_n = nucleus_table->getExBinN(target.c_str(),i);
      const int numofexbin_p = nucleus_table->getExBinP(target.c_str(),i);
      const int numofexbin_a = nucleus_table->getExBinA(target.c_str(),i);
      const int numofexbin_d = nucleus_table->getExBinD(target.c_str(),i);
      const int numofexbin_t = nucleus_table->getExBinT(target.c_str(),i);
      const int numofexbin_h = nucleus_table->getExBinH(target.c_str(),i);
      int index_g=0, index_n=0, index_p=0, index_a=0, index_d=0, index_t=0, index_h=0;
      if(total_p_g>0){
        for(int j=0; j<=numofexbin_g ;j++){
          double e_g = nucleus_table->getExG(target.c_str(),i,j);
          double p_g = nucleus_table->getPopG(target.c_str(),i,j);
          g_target_br_ex_g[i]->SetPoint(index_g,e_g,p_g/total_p_g);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_g << "   "
                  << std::setw(5) << e_g << "   " << std::setw(9) << p_g << std::endl;
          index_g++;
        }
      }else{
          g_target_br_ex_g[i]->SetPoint(index_g,0,0);
          index_g++;
      }
      if(total_p_n>0){
        for(int j=0; j<=numofexbin_n ;j++){
          double e_n = nucleus_table->getExN(target.c_str(),i,j);
          double p_n = nucleus_table->getPopN(target.c_str(),i,j);
          g_target_br_ex_n[i]->SetPoint(index_n,e_n,p_n/total_p_n);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_n << "   "
                  << std::setw(5) << e_n << "   " << std::setw(9) << p_n << std::endl;
          index_n++;
        }
      }else{
          g_target_br_ex_n[i]->SetPoint(index_n,0,0);
          index_n++;
      }
      if(total_p_p>0){
        for(int j=0; j<=numofexbin_p ;j++){
          double e_p = nucleus_table->getExP(target.c_str(),i,j);
          double p_p = nucleus_table->getPopP(target.c_str(),i,j);
          g_target_br_ex_p[i]->SetPoint(index_p,e_p,p_p/total_p_p);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_p << "   "
                  << std::setw(5) << e_p << "   " << std::setw(9) << p_p << std::endl;
          index_p++;
        }
      }else{
          g_target_br_ex_p[i]->SetPoint(index_p,0,0);
          index_p++;
      }
      if(total_p_a>0){
        for(int j=0; j<=numofexbin_a ;j++){
          double e_a = nucleus_table->getExA(target.c_str(),i,j);
          double p_a = nucleus_table->getPopA(target.c_str(),i,j);
          g_target_br_ex_a[i]->SetPoint(index_a,e_a,p_a/total_p_a);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_a << "   "
                  << std::setw(5) << e_a << "   " << std::setw(9) << p_a << std::endl;
          index_a++;
        }
      }else{
          g_target_br_ex_a[i]->SetPoint(index_a,0,0);
          index_a++;
      }
      if(total_p_d>0){
        for(int j=0; j<=numofexbin_d ;j++){
          double e_d = nucleus_table->getExD(target.c_str(),i,j);
          double p_d = nucleus_table->getPopD(target.c_str(),i,j);
          g_target_br_ex_d[i]->SetPoint(index_d,e_d,p_d/total_p_d);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_d << "   "
                  << std::setw(5) << e_d << "   " << std::setw(9) << p_d << std::endl;
          index_d++;
        }
      }else{
          g_target_br_ex_d[i]->SetPoint(index_d,0,0);
          index_d++;
      }
      if(total_p_t>0){
        for(int j=0; j<=numofexbin_t ;j++){
          double e_t = nucleus_table->getExT(target.c_str(),i,j);
          double p_t = nucleus_table->getPopT(target.c_str(),i,j);
          g_target_br_ex_t[i]->SetPoint(index_t,e_t,p_t/total_p_t);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_t << "   "
                  << std::setw(5) << e_t << "   " << std::setw(9) << p_t << std::endl;
          index_t++;
        }
      }else{
          g_target_br_ex_t[i]->SetPoint(index_t,0,0);
          index_t++;
      }
      if(total_p_h>0){
        for(int j=0; j<=numofexbin_h ;j++){
          double e_h = nucleus_table->getExH(target.c_str(),i,j);
          double p_h = nucleus_table->getPopH(target.c_str(),i,j);
          g_target_br_ex_h[i]->SetPoint(index_h,e_h,p_h/total_p_h);
          std::cout << std::setw(5) << e << "   " << std::setw(9) << total_p_h << "   "
                  << std::setw(5) << e_h << "   " << std::setw(9) << p_h << std::endl;
          index_h++;
        }
      }else{
          g_target_br_ex_h[i]->SetPoint(index_h,0,0);
          index_h++;
      }
      
      //
      //std::cout << std::setw(5) << e << "   " << std::setw(9) << p << "   "
      //          << std::setw(9) << total_p_g  << "   "
      //          << std::setw(9) << total_p_n  << "   "
      //          << std::setw(9) << total_p_p  << "   "
      //          << std::setw(9) << total_p_a  << "   " 
      //          << std::endl;
      index++;
      if(flag_below_S==1) continue;


      TCanvas* c_target_br_ex = new TCanvas("c_target_br_ex","",0,0,3200,1200);
      c_target_br_ex->Divide(4,2);
      //
      os.str("");
      os << target.c_str() << ", (bin,Ex) = (" << i << ", " << e << ")";
      TText* text_target_br_ex = new TText(30,0.9,os.str().c_str());
      text_target_br_ex->SetTextSize(0.05);
      //
      c_target_br_ex->cd(1);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_g  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_g->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_g->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(g) = " << std::setw(5) << total_p_g/p;
      waku_target_br_ex_g->SetTitle(os.str().c_str());
      g_target_br_ex_g[i]->SetLineWidth(2);
      g_target_br_ex_g[i]->SetMarkerStyle(7);
      g_target_br_ex_g[i]->SetMarkerColor(kViolet);
      g_target_br_ex_g[i]->SetLineColor(kViolet);
      g_target_br_ex_g[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      c_target_br_ex->cd(2);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_n  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_n->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_n->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(n) = " << std::setw(5) << total_p_n/p;
      waku_target_br_ex_n->SetTitle(os.str().c_str());
      g_target_br_ex_n[i]->SetLineWidth(2);
      g_target_br_ex_n[i]->SetMarkerStyle(7);
      g_target_br_ex_n[i]->SetMarkerColor(kBlue);
      g_target_br_ex_n[i]->SetLineColor(kBlue);
      g_target_br_ex_n[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      c_target_br_ex->cd(3);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_p  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_p->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_p->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(p) = " << std::setw(5) << total_p_p/p;
      waku_target_br_ex_p->SetTitle(os.str().c_str());
      g_target_br_ex_p[i]->SetLineWidth(2);
      g_target_br_ex_p[i]->SetMarkerStyle(7);
      g_target_br_ex_p[i]->SetMarkerColor(kRed);
      g_target_br_ex_p[i]->SetLineColor(kRed);
      g_target_br_ex_p[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      c_target_br_ex->cd(4);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_a  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_a->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_a->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(a) = " << std::setw(5) << total_p_a/p;
      waku_target_br_ex_a->SetTitle(os.str().c_str());
      g_target_br_ex_a[i]->SetLineWidth(2);
      g_target_br_ex_a[i]->SetMarkerStyle(7);
      g_target_br_ex_a[i]->SetMarkerColor(kMagenta);
      g_target_br_ex_a[i]->SetLineColor(kMagenta);
      g_target_br_ex_a[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      c_target_br_ex->cd(5);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_d  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_d->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_d->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(d) = " << std::setw(5) << total_p_d/p;
      waku_target_br_ex_d->SetTitle(os.str().c_str());
      g_target_br_ex_d[i]->SetLineWidth(2);
      g_target_br_ex_d[i]->SetMarkerStyle(7);
      g_target_br_ex_d[i]->SetMarkerColor(kGray+1);
      g_target_br_ex_d[i]->SetLineColor(kGray+1);
      g_target_br_ex_d[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      c_target_br_ex->cd(6);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_t  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_t->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_t->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(t) = " << std::setw(5) << total_p_t/p;
      waku_target_br_ex_t->SetTitle(os.str().c_str());
      g_target_br_ex_t[i]->SetLineWidth(2);
      g_target_br_ex_t[i]->SetMarkerStyle(7);
      g_target_br_ex_t[i]->SetMarkerColor(kCyan+1);
      g_target_br_ex_t[i]->SetLineColor(kCyan+1);
      g_target_br_ex_t[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      c_target_br_ex->cd(7);
      gPad->SetGrid();
      TH1F* waku_target_br_ex_h  = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
      waku_target_br_ex_h->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_br_ex_h->GetYaxis()->SetTitle("Branching ratio");
      os.str("");
      os << "Br(h) = " << std::setw(5) << total_p_h/p;
      waku_target_br_ex_h->SetTitle(os.str().c_str());
      g_target_br_ex_h[i]->SetLineWidth(2);
      g_target_br_ex_h[i]->SetMarkerStyle(7);
      g_target_br_ex_h[i]->SetMarkerColor(kGreen+1);
      g_target_br_ex_h[i]->SetLineColor(kGreen+1);
      g_target_br_ex_h[i]->Draw("PLsame");
      //text_target_br_ex->Draw("same");
      //
      os.str("");
      os << "fig/fig_target_br_Ex_" << i << "_" << e << "_" << target.c_str() << ".pdf";
      c_target_br_ex->Print(os.str().c_str());
      //
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_g_" << i;
      g_target_br_ex_g[i]->SetName(os.str().c_str());
      g_target_br_ex_g[i]->Write();
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_n_" << i;
      g_target_br_ex_n[i]->SetName(os.str().c_str());
      g_target_br_ex_n[i]->Write();
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_p_" << i;
      g_target_br_ex_p[i]->SetName(os.str().c_str());
      g_target_br_ex_p[i]->Write();
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_a_" << i;
      g_target_br_ex_a[i]->SetName(os.str().c_str());
      g_target_br_ex_a[i]->Write();
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_d_" << i;
      g_target_br_ex_d[i]->SetName(os.str().c_str());
      g_target_br_ex_d[i]->Write();
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_t_" << i;
      g_target_br_ex_t[i]->SetName(os.str().c_str());
      g_target_br_ex_t[i]->Write();
      os.str("");
      os << "g_" << target.c_str() << "_br_ex_h_" << i;
      g_target_br_ex_h[i]->SetName(os.str().c_str());
      g_target_br_ex_h[i]->Write();
      //
      c_target_br_ex->Update();
      c_target_br_ex->Clear();
      delete c_target_br_ex;
    }
    std::cout << total_p << std::endl;


    TCanvas* c_target_pop = new TCanvas("c_target_pop","",0,0,800,600);
    TH1F* waku_target_pop = c_target_pop->DrawFrame(0,0,70,max_p*1.1);
    os.str("");
    os << "Population (normarized), " << target.c_str();
    waku_target_pop->SetTitle(os.str().c_str());
    waku_target_pop->GetXaxis()->SetTitle("Excitation energy [MeV]");
    waku_target_pop->GetYaxis()->SetTitle("Population");
    gPad->SetGrid();
    g_target_pop->SetLineWidth(2);
    g_target_pop->SetMarkerStyle(7);
    g_target_pop->Draw("PLsame");
    os.str("");
    os << "fig/fig_target_pop_" << target.c_str() << ".pdf";
    c_target_pop->Print(os.str().c_str());
    c_target_pop->Update();
    c_target_pop->Clear();
    delete c_target_pop;
    //
    os.str("");
    os << "g_" << target.c_str() << "_pop";
    g_target_pop->SetName(os.str().c_str());
    g_target_pop->Write();



    TCanvas* c_target_br = new TCanvas("c_target_br","",0,0,800,600);
    TH1F* waku_target_br = c_target_br->DrawFrame(0,0,max_Ex_plot,1.05);
    gPad->SetGrid();
    os.str("");
    os << "Branching ratio, " << target.c_str();
    waku_target_br->SetTitle(os.str().c_str());
    waku_target_br->GetXaxis()->SetTitle("Excitation energy [MeV]");
    waku_target_br->GetYaxis()->SetTitle("Branching ratio");
    waku_target_br->GetXaxis()->SetTitleSize(0.045);
    waku_target_br->GetYaxis()->SetTitleSize(0.045);
    //
    g_target_br_g->SetLineWidth(2);
    g_target_br_g->SetMarkerStyle(7);
    g_target_br_g->SetLineColor(kViolet);
    g_target_br_g->SetMarkerColor(kViolet);
    g_target_br_g->Draw("PLsame");
    //
    g_target_br_n->SetLineWidth(2);
    g_target_br_n->SetMarkerStyle(7);
    g_target_br_n->SetLineColor(kBlue);
    g_target_br_n->SetMarkerColor(kBlue);
    g_target_br_n->Draw("PLsame");
    //
    g_target_br_p->SetLineWidth(2);
    g_target_br_p->SetMarkerStyle(7);
    g_target_br_p->SetLineColor(kRed);
    g_target_br_p->SetMarkerColor(kRed);
    g_target_br_p->Draw("PLsame");
    //
    g_target_br_a->SetLineWidth(2);
    g_target_br_a->SetMarkerStyle(7);
    g_target_br_a->SetLineColor(kMagenta);
    g_target_br_a->SetMarkerColor(kMagenta);
    g_target_br_a->Draw("PLsame");
    //
    g_target_br_d->SetLineWidth(2);
    g_target_br_d->SetMarkerStyle(7);
    g_target_br_d->SetLineColor(kGray+1);
    g_target_br_d->SetMarkerColor(kGray+1);
    g_target_br_d->Draw("PLsame");
    //
    g_target_br_t->SetLineWidth(2);
    g_target_br_t->SetMarkerStyle(7);
    g_target_br_t->SetLineColor(kCyan+1);
    g_target_br_t->SetMarkerColor(kCyan+1);
    g_target_br_t->Draw("PLsame");
    //
    g_target_br_h->SetLineWidth(2);
    g_target_br_h->SetMarkerStyle(7);
    g_target_br_h->SetLineColor(kGreen+1);
    g_target_br_h->SetMarkerColor(kGreen+1);
    g_target_br_h->Draw("PLsame");
    //
    os.str("");
    os << "fig/fig_target_br_" << target.c_str() << ".pdf";
    c_target_br->Print(os.str().c_str());
    c_target_br->Update();
    c_target_br->Clear();
    delete c_target_br;
    //
    os.str("");
    os << "g_" << target.c_str() << "_br_g";
    g_target_br_g->SetName(os.str().c_str());
    g_target_br_g->Write();
    os.str("");
    os << "g_" << target.c_str() << "_br_n";
    g_target_br_n->SetName(os.str().c_str());
    g_target_br_n->Write();
    os.str("");
    os << "g_" << target.c_str() << "_br_p";
    g_target_br_p->SetName(os.str().c_str());
    g_target_br_p->Write();
    os.str("");
    os << "g_" << target.c_str() << "_br_a";
    g_target_br_a->SetName(os.str().c_str());
    g_target_br_a->Write();
    os.str("");
    os << "g_" << target.c_str() << "_br_d";
    g_target_br_d->SetName(os.str().c_str());
    g_target_br_d->Write();
    os.str("");
    os << "g_" << target.c_str() << "_br_t";
    g_target_br_t->SetName(os.str().c_str());
    g_target_br_t->Write();
    os.str("");
    os << "g_" << target.c_str() << "_br_h";
    g_target_br_h->SetName(os.str().c_str());
    g_target_br_h->Write();
    rootf->Close();
    delete rootf; 

  // output, target
    std::ofstream ofs;
    //
    os.str("");
    os << "txt/Br_" << target.c_str() << "_summary.txt";
    ofs.open(os.str().c_str());
    //
    Sg=-1, Sn=-1, Sp=-1, Sd=-1, St=-1, Sh=-1, Sa=-1;
    os.str("");
    os << "/userdata/work/seisho/TALYS/work-dir_talys_1.95/separation_energy/separation_energy_" << target.c_str() <<  ".txt";
    ifs.open(os.str().c_str());
    if(!ifs.is_open()) return 1;
    else std::cout << os.str().c_str() << std::endl;
    while(ifs.getline(buf,sizeof(buf))){
      //std::istringstream(buf) >> Sg >> Sn >> Sp >> Sd >> St >> Sh >> Sa;
      if(!( std::istringstream(buf) >> Sg >> Sn >> Sp >> Sd >> St >> Sh >> Sa)){
        cerr << "Something wrong in separation E table " << endl;
        return 1;
      }
    }
    ifs.close();
    //
    for(int i=0;i<=Ex_max;i+=Ex_bin_width){ //keV
      bool flag_g=0, flag_n=0, flag_p=0, flag_a=0, flag_d=0, flag_t=0, flag_h=0;
      double Br_g = g_target_br_g->Eval(i*1e-3); // keV2MeV
      double Br_n = g_target_br_n->Eval(i*1e-3); // keV2MeV
      double Br_p = g_target_br_p->Eval(i*1e-3); // keV2MeV
      double Br_a = g_target_br_a->Eval(i*1e-3); // keV2MeV
      double Br_d = g_target_br_d->Eval(i*1e-3); // keV2MeV
      double Br_t = g_target_br_t->Eval(i*1e-3); // keV2MeV
      double Br_h = g_target_br_h->Eval(i*1e-3); // keV2MeV
      if(Br_g>0 && i*1e-3>Sg) flag_g=1;
      if(Br_n>0 && i*1e-3>Sn) flag_n=1;
      if(Br_p>0 && i*1e-3>Sp) flag_p=1;
      if(Br_a>0 && i*1e-3>Sa) flag_a=1;
      if(Br_d>0 && i*1e-3>Sd) flag_d=1;
      if(Br_t>0 && i*1e-3>St) flag_t=1;
      if(Br_h>0 && i*1e-3>Sh) flag_h=1;
      if(flag_n==0 && flag_p==0 && flag_a==0
          && flag_d==0 && flag_t==0 && flag_h==0) continue;
      ofs << "P     "  << std::setw(15) << resetiosflags(std::ios_base::floatfield) << i << "  - 0" << std::endl; 
      if(flag_g==1) ofs << "                    IT            0     " << Br_g << std::endl;
      if(flag_n==1) ofs << "                    Neutron       0     " << Br_n << std::endl;
      if(flag_p==1) ofs << "                    Proton        0     " << Br_p << std::endl;
      if(flag_a==1) ofs << "                    Alpha         0     " << Br_a << std::endl;
      if(flag_d==1) ofs << "                    Deuteron      0     " << Br_d << std::endl;
      if(flag_t==1) ofs << "                    Triton        0     " << Br_t << std::endl;
      if(flag_h==1) ofs << "                    He3           0     " << Br_h << std::endl;
      // get nearest ex  (common for decay type)
      double X, X_tmp;
      int bin_tmp;
      for(int j=0;j<g_target_br_g->GetN();j++){
        double Y;
        g_target_br_g->GetPoint(j,X,Y);
        //if(X>i*1e-3) break;
        if(X>i*1e-3){
          X=X_tmp;
          break;
        }
        X_tmp=X;
      }
      for(int j=0;j<nucleus_table->getExBin(target.c_str());j++){
        if(X==nucleus_table->getEx(target.c_str(),j)){
          bin_tmp=j;
          break;
        }
      }
      if(flag_n==1){
        // get br for each level of daughter 
        double norm=0.;
        for(int j=0;j<g_target_br_ex_n[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_n[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_n) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_n[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sn - Xd)*1e3; // keV
          if(Q>0) norm += Yd;
        }
        for(int j=0;j<g_target_br_ex_n[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_n[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_n) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_n[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sn - Xd)*1e3; // keV
          double ratio = Br_n*100*Yd/norm; // already normarized
          if(Q>0){
            ofs << "                    Neutron       " << std::setw(9) << Xd*1e3 << "   -   " << std::setw(15) << ratio << "   "
                 << std::setw(15) << Q << std::endl;
          }
        }
      }
      if(flag_p==1){
        // get br for each level of daughter 
        double norm=0.;
        for(int j=0;j<g_target_br_ex_p[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_p[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_p) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_p[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sp - Xd)*1e3; // keV
          if(Q>0) norm += Yd;
        }
        for(int j=0;j<g_target_br_ex_p[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_p[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_p) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_p[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sp - Xd)*1e3; // keV
          double ratio = Br_p*100*Yd/norm; // already normarized
          if(Q>0){
            ofs << "                     Proton       " << std::setw(9) << Xd*1e3 << "   -   " << std::setw(15) << ratio << "   "
                 << std::setw(15) << Q << std::endl;
          }
        }
      }
      if(flag_a==1){
        // get br for each level of daughter 
        double norm=0.;
        for(int j=0;j<g_target_br_ex_a[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_a[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_a) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_a[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sa - Xd)*1e3; // keV
          if(Q>0) norm += Yd;
        }
        for(int j=0;j<g_target_br_ex_a[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_a[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_a) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_a[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sa - Xd)*1e3; // keV
          double ratio = Br_a*100*Yd/norm; // already normarized
          if(Q>0){
            ofs << "                     Alpha        " << std::setw(9) << Xd*1e3 << "   -   " << std::setw(15) << ratio << "   "
                 << std::setw(15) << Q << std::endl;
          }
        }
      }
      if(flag_d==1){
        // get br for each level of daughter 
        double norm=0.;
        for(int j=0;j<g_target_br_ex_d[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_d[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_d) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_d[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sd - Xd)*1e3; // keV
          if(Q>0) norm += Yd;
        }
        for(int j=0;j<g_target_br_ex_d[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_d[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_d) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_d[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sd - Xd)*1e3; // keV
          double ratio = Br_d*100*Yd/norm; // already normarized
          if(Q>0){
            ofs << "                     Deuteron     " << std::setw(9) << Xd*1e3 << "   -   " << std::setw(15) << ratio << "   "
                 << std::setw(15) << Q << std::endl;
          }
        }
      }
      if(flag_t==1){
        // get br for each level of daughter 
        double norm=0.;
        for(int j=0;j<g_target_br_ex_t[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_t[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_t) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_t[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - St - Xd)*1e3; // keV
          if(Q>0) norm += Yd;
        }
        for(int j=0;j<g_target_br_ex_t[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_t[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_t) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_t[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - St - Xd)*1e3; // keV
          double ratio = Br_t*100*Yd/norm; // already normarized
          if(Q>0){
            ofs << "                     Triton       " << std::setw(9) << Xd*1e3 << "   -   " << std::setw(15) << ratio << "   "
                 << std::setw(15) << Q << std::endl;
          }
        }
      }
      if(flag_h==1){
        // get br for each level of daughter 
        double norm=0.;
        for(int j=0;j<g_target_br_ex_h[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_h[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_h) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_h[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sh - Xd)*1e3; // keV
          if(Q>0) norm += Yd;
        }
        for(int j=0;j<g_target_br_ex_h[bin_tmp]->GetN();j++){
          double Xd, Yd;
          g_target_br_ex_h[bin_tmp]->GetPoint(j,Xd,Yd);
          if(j>maxlevelsbin_h) Xd = ((int)(Xd*1.0e3)/Ex_bin_width)*Ex_bin_width/1.0e3;
          if(g_target_br_ex_h[bin_tmp]->GetN()==1 && Yd==0) Yd=1.0;
          double Q = (i*1e-3 - Sh - Xd)*1e3; // keV
          double ratio = Br_h*100*Yd/norm; // already normarized
          if(Q>0){
            ofs << "                     He3          " << std::setw(9) << Xd*1e3 << "   -   " << std::setw(15) << ratio << "   "
                 << std::setw(15) << Q << std::endl;
          }
        }
      }
    }//  end of ex lopo
    ofs.close();
  }

  return 0;
}
