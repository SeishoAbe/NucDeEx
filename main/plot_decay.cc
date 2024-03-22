#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "NucDeExConsts.hh"
#include "NucDeExUtils.hh"
#include "NucDeExNucleusTable.hh"
#include "ReadTALYS.hh"

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TTree.h>


const std::string decay_name[NucDeEx::num_particle] // for G4
  = {"IT","Neutron","Proton",
     "Deuteron","Triton","He3","Alpha"};

int main(int argc, char* argv[]){
  if(argc!=5){
    std::cerr << argv[0] << " [Target nucleus] [ldmodel] [parity&optmodall] [flag_jpi]" << std::endl;
    return 0;
  }
  // --- FIXME  --- //
  const float max_Ex_plot=50; // <- plot
  // ---------------//

  const int ldmodel = atoi(argv[2]);
  const bool parity_optmodall = (bool)atoi(argv[3]);
  const bool flag_jpi = (bool) atoi(argv[4]);
  
  std::ostringstream os;
  NucDeExNucleusTable* nucleus_table = new NucDeExNucleusTable();
  if(!nucleus_table->ReadTables()){
    std::cerr << "something wrong" << std::endl;
    return 1;
  }
  NucDeExNucleus* nuc = nucleus_table->GetNucleusPtr(argv[1]);
  const int Zt=nuc->Z;
  const int Nt=nuc->N;
  const int At=Zt+Nt;

  os.str("");
  os << "output/";
  if(flag_jpi){
    if(At==11) os << "12C";
    else if(At==15) os << "16O";
    else abort();
  }
  os << "/output_" << argv[1] << "_ldmodel" << ldmodel;
  if(parity_optmodall) os << "_parity_optmodall";
  NucDeExUtils::SetVerbose(1);
  ReadTALYS* read_talys = new ReadTALYS(os.str().c_str(), nucleus_table);
  if(!read_talys->Read()){
    std::cerr << "something wrong happend..." << std::endl;
    return 0;
  }

  // Get general information
  //
  const int num_of_nuc = nucleus_table->GetNumofNuc();
  int num_nuc_data=0;
  for(int inuc=0;inuc<num_of_nuc;inuc++){
    if(nucleus_table->GetNucleusPtr(inuc)->flag_data==1)  num_nuc_data++;
  }
  std::cout << "# nuclei that has population data = " << num_nuc_data << std::endl;
  

  // prepare output root file
  os.str("");
  os << "output/";
  if(flag_jpi){
    if(At==11) os << "12C";
    else if(At==15) os << "16O";
    else abort();
  }
  os << "/Br_" << argv[1] << "_ldmodel" << ldmodel;
  if(parity_optmodall) os << "_parity_optmodall";
  os << ".root";
  TFile* rootf = new TFile(os.str().c_str(),"RECREATE");
  TTree* tree = new TTree("tree","tree");
  std::string name;
  int Z,N;
  int Ex_bin[NucDeEx::num_particle];
  float Pop[NucDeEx::parity][NucDeEx::bins];
  float Ex[NucDeEx::num_particle][NucDeEx::bins];
  float Br[NucDeEx::num_particle][NucDeEx::bins];
  int REx_bin[NucDeEx::num_particle][NucDeEx::bins];
  float REx[NucDeEx::num_particle][NucDeEx::bins][NucDeEx::bins], RBr[NucDeEx::num_particle][NucDeEx::bins][NucDeEx::bins];
  tree->Branch("name",&name);
  tree->Branch("Z",&Z,"Z/I");
  tree->Branch("N",&N,"N/I");
  os.str("");
  os << "Pop[" << NucDeEx::parity << "][" << NucDeEx::bins << "]/F";
  tree->Branch("Pop",&Pop,os.str().c_str());
  os.str("");
  os << "Ex_bin[" << NucDeEx::num_particle << "]/I";
  tree->Branch("Ex_bin",&Ex_bin,os.str().c_str());
  os.str("");
  os << "Ex[" << NucDeEx::num_particle << "][" << NucDeEx::bins << "]/F";
  tree->Branch("Ex",&Ex,os.str().c_str());
  os.str("");
  os << "Br[" << NucDeEx::num_particle << "][" << NucDeEx::bins << "]/F";
  tree->Branch("Br",&Br,os.str().c_str());
  //
  os.str("");
  os << "REx_bin[" << NucDeEx::num_particle << "][" << NucDeEx::bins << "]/I";
  tree->Branch("REx_bin",&REx_bin,os.str().c_str());
  os.str("");
  os << "REx[" << NucDeEx::num_particle << "][" << NucDeEx::bins << "][" << NucDeEx::bins << "]/F";
  tree->Branch("REx",&REx,os.str().c_str());
  os.str("");
  os << "RBr[" << NucDeEx::num_particle << "][" << NucDeEx::bins << "][" << NucDeEx::bins << "]/F";
  tree->Branch("RBr",&RBr,os.str().c_str());

  // start Nucleus loop in NucleusTable
  for(int inuc=0;inuc<num_of_nuc;inuc++){
    if(!nucleus_table->GetNucleusPtr(inuc)->flag_data) continue; // no data -> continue

    NucDeExNucleus* nuc_target = nucleus_table->GetNucleusPtr(inuc); // get pointer of Nulceus

    // check info in Nucleus 
    if(! (nuc_target->CheckTotalPop() && nuc_target->CheckPop()) ){
      std::cerr << "ERROR: There was a bug in population info: " << nuc_target->name << std::endl;
      return 0;
    }else{
      std::cout << "Population for " << nuc_target->name << " looks OK" << std::endl;
    }
    if(! (nuc_target->CheckEx())) {
      std::cerr << "ERROR: There was unexpected behaviour in Ex" << std::endl;
      return 0;
    }

    // Init array
    for(int p=0;p<NucDeEx::num_particle;p++){
      Ex_bin[p]=0;
      for(int b=0;b<NucDeEx::bins;b++){
        Ex[p][b]=Br[p][b]=0;
        REx_bin[p][b]=0;
        for(int b1=0;b1<NucDeEx::bins;b1++){
          REx[p][b][b1]=RBr[p][b][b1]=0;
        }
      }
    }

    // for ttree
    name = (std::string)nuc_target->name;
    Z = nuc_target->Z;
    N = nuc_target->N;
    int bin_target = nuc_target->Ex_bin[0];

    int index[NucDeEx::parity]={0};
    int  index_br[NucDeEx::num_particle]={0};
    int index_br_ex[NucDeEx::num_particle][NucDeEx::bins]={0};
    float max_total_pop[NucDeEx::parity]={0};
    for(int par=0;par<NucDeEx::parity;par++){
      index[par]=0;
      max_total_pop[par]=0;
    }
    for(int p=0;p<NucDeEx::num_particle;p++){
      index_br[p] = 0;
      for(int i=0;i<NucDeEx::bins;i++){ // target bin loop
        index_br_ex[p][i]=0;
      }
    }

    for(int i=0;i<bin_target;i++){ // target bin loop
      float Ex_target = nuc_target->Ex[0][i];
      float pop_target_sum = nuc_target->GetPopParitySum(i);

      for(int par=0;par<NucDeEx::parity;par++){
        float pop_target_norm = nuc_target->pop[par][i]/nuc_target->total_pop[par];
        // fill normalized population
        //g_target_pop[par]->SetPoint(index[par],Ex_target,pop_target_norm);
        Pop[par][index[par]] = pop_target_norm;
        if(max_total_pop[par]<pop_target_norm) max_total_pop[par] = pop_target_norm;
        index[par]++;
      }

      // calculate br
      for(int p=0;p<NucDeEx::num_particle;p++){ // particle loop
        if(!nuc_target->flag_decay_data[i] && i<=20){ // no decay data
          if(p==0){
            //g_target_br[p]->SetPoint(index_br[p],Ex_target,1);
            Ex[p][index_br[p]] = Ex_target;
            Br[p][index_br[p]] = 1;
            // Fill level info only for gamma
            for(int li=0;li<=i;li++){
              //g_target_br_ex[p][i]->SetPoint(index_br_ex[p][i],  
              //                               nuc_target->level_Ex[li],
              //                               nuc_target->level_br[i][li]*1e-2); // % to frac
              REx[p][index_br[p]][index_br_ex[p][index_br[p]]] = nuc_target->level_Ex[li];
              RBr[p][index_br[p]][index_br_ex[p][index_br[p]]] = nuc_target->level_br[i][li]*1e-2;
              index_br_ex[p][index_br[p]]++;
            }
            index_br[p]++;
          }else{
            //g_target_br[p]->SetPoint(index_br[p],Ex_target,0);
            Ex[p][index_br[p]] = Ex_target;
            Br[p][index_br[p]] = 0;
            index_br[p]++;
          }
        }else if(pop_target_sum>0 && nuc_target->flag_decay_data[i]){ // have decay data
          float population = nuc_target->GetPopDaughterBinSum(p,i);//  (particle, exbin)
          float br = population/pop_target_sum;
          if(br<=0) br=0;
          if(br>=1) br=1;
          //g_target_br[p]->SetPoint(index_br[p],Ex_target,br);
          Ex[p][index_br[p]] = Ex_target;
          Br[p][index_br[p]] = br;
          for(int j=0;j<nuc_target->Ex_bin_p[p][i] && population>0;j++){ // daughter ex bin loop
            float br_ex = nuc_target->pop_p[p][i][j]/population;
            if(br_ex<=0) br_ex=0;
            if(br_ex>=1) br_ex=1;
            //g_target_br_ex[p][i]->SetPoint(index_br_ex[p][i],
            //                               nuc_target->Ex_p[p][i][j],br_ex);
            REx[p][index_br[p]][index_br_ex[p][index_br[p]]] = nuc_target->Ex_p[p][i][j];
            RBr[p][index_br[p]][index_br_ex[p][index_br[p]]] = br_ex;
            index_br_ex[p][index_br[p]]++;
          }
          index_br[p]++;
        }else{
          std::cerr << "Warning: No decay data, but have population " << std::endl;
          std::cerr << "pop_targe_sum = " << pop_target_sum << std::endl;
          std::cerr << "decay_data = " << nuc_target->flag_decay_data[i] << std::endl;
          std::cerr << "bin = " << i << std::endl;
          std::cerr << "particle = " << p << std::endl;
        }
      } // end of p (particle) loop
    }// end of i (target bin) loop

    for(int p=0;p<NucDeEx::num_particle;p++){ // particle loop
      Ex_bin[p]       = index_br[p];
      for(int i=0;i<index_br[p];i++){ // target bin loop
        REx_bin[p][i]  = index_br_ex[p][i];
      }
    }

    // --- Fill Tree --- //
    tree->Fill();

    // --- Prepare TGraph --- //
    TGraph* g_target_pop[NucDeEx::parity];
    TGraph* g_target_br[NucDeEx::num_particle];
    TGraph* g_target_br_ex[NucDeEx::num_particle][bin_target];
    for(int par=0;par<NucDeEx::parity;par++){
      g_target_pop[par] = new TGraph(index[par],Ex[0],Pop[par]);
      os.str("");
      os << "g_" << name.c_str() << "_pop_" << par;
      g_target_pop[par]->SetName(os.str().c_str());
    }
    for(int p=0;p<NucDeEx::num_particle;p++){
      g_target_br[p] = new TGraph(Ex_bin[p],Ex[p],Br[p]);
      os.str("");
      os << "g_" << name.c_str() << "_br_" << p;
      g_target_br[p]->SetName(os.str().c_str());
      g_target_br[p]->SetMarkerStyle(7);
      g_target_br[p]->SetMarkerColor(NucDeEx::color_root[p]);
      g_target_br[p]->SetLineWidth(2);
      g_target_br[p]->SetLineColor(NucDeEx::color_root[p]);
    }
    for(int p=0;p<NucDeEx::num_particle;p++){
      for(int i=0;i<bin_target;i++){ // target bin loop
        g_target_br_ex[p][i] = new TGraph(REx_bin[p][i],REx[p][i],RBr[p][i]);
        os.str("");
        os << "g_" << name << "_br_ex_" << p << "_" << i;
        g_target_br_ex[p][i]->SetName(os.str().c_str());
        g_target_br_ex[p][i]->SetMarkerStyle(7);
        g_target_br_ex[p][i]->SetMarkerColor(NucDeEx::color_root[p]);
        g_target_br_ex[p][i]->SetLineWidth(2);
        g_target_br_ex[p][i]->SetLineColor(NucDeEx::color_root[p]);
      }
    }


    // --- Draw --- // 
    gStyle->SetTextSize(0.08);
    gStyle->SetTitleSize(0.045);
    gStyle->SetTitleXSize(0.045);
    gStyle->SetTitleYSize(0.045);
    gStyle->SetTitleYOffset(0.95);
    
    // --- Canvas --- //
    TCanvas* c_target_pop = new TCanvas("c_target_pop","",0,0,1200,600);
    c_target_pop->Divide(2);
    for(int par=0;par<NucDeEx::parity;par++){
      c_target_pop->cd(par+1);
      TH1F* waku_target_pop = gPad->DrawFrame(0,0,max_Ex_plot,max_total_pop[par]*1.1);
      os.str("");
      os << "Population (normarized), " << name.c_str();
      waku_target_pop->SetTitle(os.str().c_str());
      waku_target_pop->GetXaxis()->SetTitle("Excitation energy [MeV]");
      waku_target_pop->GetYaxis()->SetTitle("Population (A.U.)");
      gPad->SetGrid();
      g_target_pop[par]->SetLineWidth(2);
      g_target_pop[par]->SetMarkerStyle(7);
      g_target_pop[par]->Draw("PLsame");
    }
    gPad->RedrawAxis();
    os.str("");
    os << "fig/";
    if(flag_jpi){
      if(At==11) os << "12C/";
      else if(At==15) os << "16O/";
      else abort();
    }
    os << argv[1] << "_ldmodel" << ldmodel;
    if(parity_optmodall) os << "_parity_optmodall";
    os << "/fig_" << name.c_str() << "_pop.pdf";
    c_target_pop->Print(os.str().c_str());
    c_target_pop->Update();
    c_target_pop->Clear();
    delete c_target_pop;



    TCanvas* c_target_br = new TCanvas("c_target_br","",0,0,800,600);
    TH1F* waku_target_br = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
    os.str("");
    os << "Branching ratio: target: " << name.c_str();
    waku_target_br->SetTitle(os.str().c_str());
    waku_target_br->GetXaxis()->SetTitle("Excitation energy [MeV]");
    waku_target_br->GetYaxis()->SetTitle("Branching ratio");
    gPad->SetGrid();
    for(int p=0;p<NucDeEx::num_particle;p++){
      g_target_br[p]->Draw("PLsame");
    }
    TLegend* leg_target_br = new TLegend(0.4,0.7,0.9,0.9);
    leg_target_br->SetNColumns(3);
    for(int p=0;p<NucDeEx::num_particle;p++){
      leg_target_br->AddEntry(g_target_br[p],NucDeEx::particle_name[p].c_str(),"PL");
    }
    leg_target_br->Draw("same");
    gPad->RedrawAxis();
    //
    os.str("");
    os << "fig/";
    if(flag_jpi){
      if(At==11) os << "12C/";
      else if(At==15) os << "16O/";
      else abort();
    }
    os << argv[1] << "_ldmodel" << ldmodel;
    if(parity_optmodall) os << "_parity_optmodall";
    os << "/fig_" << name.c_str() << "_br.pdf";
    c_target_br->Print(os.str().c_str());
    c_target_br->Update();
    c_target_br->Clear();
    delete c_target_br;


    TCanvas* c_target_br_ex = new TCanvas("c_target_br_ex","",0,0,1200,600);
    os.str("");
    os << "fig/";
    if(flag_jpi){
      if(At==11) os << "12C/";
      else if(At==15) os << "16O/";
      else abort();
    }
    os << argv[1] << "_ldmodel" << ldmodel;
    if(parity_optmodall) os << "_parity_optmodall";
    os << "/fig_" << name.c_str() << "_br_Ex.pdf";
    std::string pdfname= os.str();
    c_target_br_ex->Print( (pdfname+(std::string)"[").c_str() );
    c_target_br_ex->Update();
    c_target_br_ex->Clear();
    for(int i=0;i<bin_target;i++){ // target bin loop
      //if(nuc_target->Ex[0][i]<nuc_target->min_S()) continue;
      float Ex_target = nuc_target->Ex[0][i];
      if(Ex_target==0) continue;
      c_target_br_ex->Divide(4,2);
      for(int p=0;p<NucDeEx::num_particle;p++){ // particle loop
        c_target_br_ex->cd(p+1);
        gPad->SetGrid();
        TH1F* waku = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
        os.str("");
        os << NucDeEx::particle_name[p];
        waku->SetTitle(os.str().c_str());
        waku->GetXaxis()->SetTitle("Excitation energy of daughter nucleus");
        waku->GetYaxis()->SetTitle("Relative branching ratio");
        if(g_target_br_ex[p][i]->GetN()==0){
          TText* text = new TText(max_Ex_plot*0.25,0.8,"NO DATA");
          text->SetTextColor(kGray+1);
          text->Draw("same");
        }else{
          g_target_br_ex[p][i]->Draw("PLsame");
          os.str("");
          //if((Ex_target<nuc_target->min_S() 
          //    || (!nuc_target->flag_decay_data[i] && i<=1)) && p==0){ 
          if((!nuc_target->flag_decay_data[i] && i<=20)){ 
            os << "BR: " << std::scientific << std::setprecision(2) << 1;
          }else{
            double br,ex;
            g_target_br[p]->GetPoint(i,ex,br);
            os << "BR: " << std::scientific << std::setprecision(2) << br;
               //<< nuc_target->GetPopDaughterBinSum(p,i)/nuc_target->GetPopParitySum(i);
          }
          TText* text = new TText(max_Ex_plot*0.25,0.8,os.str().c_str());
          text->Draw("same");
        }
      }
      c_target_br_ex->cd(8);
      os.str("");
      os << name.c_str() << ", Ex[" << i << "] = " << std::fixed << std::setprecision(2) << nuc_target->Ex[0][i] << " (MeV)";
      TText* text_Ex = new TText(0.1,0.9,os.str().c_str());
      text_Ex->Draw("same");
      TText* text_S[NucDeEx::num_particle];
      for(int p=1;p<NucDeEx::num_particle;p++){ // particle loop
        os.str("");
        os << "S_" << NucDeEx::particle_name[p].substr(0,1) << " = " << std::fixed << std::setprecision(2) << nuc_target->S[p] << " (MeV)";
        text_S[p] = new TText(0.1,0.9-0.1*p,os.str().c_str());
        if(nuc_target->S[p]>Ex_target) text_S[p]->SetTextColor(kGray+1);
        text_S[p]->Draw("same");
      }
      if((!nuc_target->flag_decay_data[i] && i<=20)){ 
        TText* text_level = new TText(0.1,0.1,"From level data");
        text_level->Draw("same");
      }
      //
      gPad->RedrawAxis();
      os.str("");
      os << pdfname.c_str();
      c_target_br_ex->Print(os.str().c_str());
      c_target_br_ex->Update();
      c_target_br_ex->Clear();
    }
    os.str("");
    os << (pdfname+(std::string)"]").c_str();
    c_target_br_ex->Print(os.str().c_str());
    c_target_br_ex->Clear();
    delete c_target_br_ex;


    // --- Save TGraph in TFile --- //
    rootf->cd();
    for(int par=0;par<NucDeEx::parity;par++){
      g_target_pop[par]->Write();
    }
    for(int p=0;p<NucDeEx::num_particle;p++){
      g_target_br[p]->Write();
      for(int i=0;i<bin_target;i++){ // target bin loop
        g_target_br_ex[p][i]->Write();
      }
    }
  }
  // end of Nucleus loop in NucleusTable

  //tree->Write();
  rootf->Close();
  delete rootf;


  delete read_talys;
  delete nucleus_table;


  return 0;
}
