#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TArrow.h>
#include <TText.h>
#include <TPaveText.h>
#include <THStack.h>

#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "NucDeExDeexcitation.hh"
#include "NucDeExEventInfo.hh"

int main(int argc, char* argv[]){
  if(argc<=5){
    std::cerr << "Input: " << argv[0] << " [Target nucleus] [ldmodel] [parity&optmodall] [flag_jpi] [verbosity] [Random Seed (optional)]" << std::endl;
    return 0;
  }
  const int ldmodel=atoi(argv[2]);
  const bool parity_optmodall=(bool)atoi(argv[3]);
  const bool flag_jpi = (bool) atoi(argv[4]);
  const int verbose = atoi(argv[5]);
  int seed=1; // default: 1
  if(argc==7) seed = atoi(argv[6]);

  std::cout << "VERBOSE = " << verbose << std::endl;
  std::cout << "SEED = " << seed << std::endl;

  // ---- FIXME --- // 
  const int numofevent=1e6; // to be generated
  const bool flag_fig=0;
  //const int numofevent=1000; // to be generated
  //const bool flag_fig=1;

  const double Ex_min =-1;
  const double Ex_max =-1;
    // negative -> not applied
  // -------------- //

  std::ostringstream os,os_remove_g;

  // --- Setup routine --- //
  // Set parameters prior to the instantiation
  NucDeEx::Utils::fVerbose=verbose; // optional (default: 0)
  NucDeEx::Random::SetSeed(seed); // optional (default: 1)
  NucDeExDeexcitation* deex = new NucDeExDeexcitation(ldmodel, parity_optmodall,2);// new v2.1~
  //NucDeExDeexcitation* deex = new NucDeExDeexcitation(ldmodel, parity_optmodall,1); // old ~v1.3
  // --------------------- //

  // Get Z and N
  NucDeExNucleus* nuc = NucDeEx::Utils::NucleusTable->GetNucleusPtr(argv[1]);
  const int Z = nuc->Z;
  const int N = nuc->N;
  const int A = Z+N;
  int Zt,Nt;
  if(Z+N==11 || Z+N==12) Zt=6,Nt=6;
  if(Z+N==15 || Z+N==16) Zt=8,Nt=8;
  
  // prepare output root file
  os.str("");  
  os << "sim_out/";
  if(flag_jpi){
    if(A==11||A==12) os << "12C/";
    else if(A==15||A==16) os << "16O/";
    else abort();
  }
  os << argv[1] << "_ldmodel" << ldmodel;
  if(parity_optmodall) os << "_parity_optmodall";
  os << ".root";
  TFile* outf = new TFile(os.str().c_str(),"RECREATE");
  TTree* tree = new TTree("tree",""); // LAB Freme, MeV
  int eventID=0, size, shell;
  double MissE, Ex, S;
  double PinitMag, PinitX,PinitY,PinitZ;
  //
  int PDG[NucDeEx::bins];
  double mass[NucDeEx::bins];
  double totalE[NucDeEx::bins],kE[NucDeEx::bins];
  double PMag[NucDeEx::bins], PX[NucDeEx::bins],PY[NucDeEx::bins],PZ[NucDeEx::bins];
  bool flag[NucDeEx::bins];
  double Ex_daughter[NucDeEx::bins];
  string decay, decay_remove_g;
  tree->Branch("eventID",&eventID,"eventID/I");
  tree->Branch("decay",&decay);
  tree->Branch("decay_remove_g",&decay_remove_g);
  tree->Branch("MissE",&MissE,"MissE/D");
  tree->Branch("S",&S,"S/D");
  tree->Branch("Ex",&Ex,"Ex/D");
  tree->Branch("shell",&shell,"shell/I");
  tree->Branch("PinitMag",&PinitMag,"PinitMag/D");
  tree->Branch("PinitX",&PinitX,"PinitX/D");
  tree->Branch("PinitY",&PinitY,"PinitY/D");
  tree->Branch("PinitZ",&PinitZ,"PinitZ/D");
  tree->Branch("size",&size,"size/I");
  tree->Branch("PDG",&PDG,"PDG[size]/I");
  tree->Branch("mass",&mass,"mass[size]/D");
  tree->Branch("totalE",&totalE,"totalE[size]/D");
  tree->Branch("kE",&kE,"kE[size]/D");
  tree->Branch("PMag",&PMag,"PMag[size]/D");
  tree->Branch("PX",&PX,"PX[size]/D");
  tree->Branch("PY",&PY,"PY[size]/D");
  tree->Branch("PZ",&PZ,"PZ[size]/D");
  tree->Branch("flag",&flag,"flag[size]/O");
  tree->Branch("Ex_daughter",&Ex_daughter,"Ex_daughter[size]/D");


  gStyle->SetTextSize(0.08);
  gStyle->SetTitleSize(0.045);
  gStyle->SetTitleXSize(0.045);
  gStyle->SetTitleYSize(0.045);
  gStyle->SetTitleYOffset(0.95);

  TCanvas* c_detail;
  string pdfname;
  if(flag_fig){
    os.str("");
    os << "fig_sim/";
    if(flag_jpi){
      if(A==11||A==12) os << "12C/";
      else if(A==15||A==16) os << "16O/";
      else abort();
    }
    os << "fig_" << argv[1] << "_ldmodel" << ldmodel;
    if(parity_optmodall) os << "_parity_optmodall";
    os << "_detail.pdf";
    pdfname = os.str();
    c_detail = new TCanvas("c_detail","",0,0,800,800);
    c_detail->Print( (pdfname + (string)"[").c_str() );
    c_detail->Update();
    c_detail->Clear();
  }
  
  //
  TVector3 Pinit;

  while(eventID<numofevent){
    // determine momentum (scalar) and missing E according to SF
    PinitMag = 0;
    Ex=25;
    // determine angle 
    double costheta = 2.*gRandom->Rndm()-1;
    double sintheta = sqrt( 1. - pow(costheta,2) );
    double phi      = 2*TMath::Pi()*gRandom->Rndm();
    Pinit.SetXYZ(PinitMag*sintheta*cos(phi),
                 PinitMag*sintheta*sin(phi),
                 PinitMag*costheta);

    // select ROI
    if(Ex_min>0 && Ex<Ex_min) continue;
    if(Ex_max>0 && Ex>Ex_max) continue;

    // --- DO SIMULATION --- //
    NucDeExEventInfo result = deex->DoDeex(Zt,Nt,Z,N,Ex,Pinit);
    /* old style ~v1.3
    NucDeExEventInfo result = deex->DoDeex(Zt,Nt,Z,N,0,Ex,Pinit); // shell level is determined from Ex (box cut)
      // 5th index
      //    1: s1/2-hole
      //    2: p3/2-hole
      //    3: p1/2-hole (nothing to do. g.s.)
    */


    // --- Scoling --- //
    shell = result.fShell;
    vector<NucDeExParticle> particle = result.ParticleVector;
    PinitX   = Pinit.X();
    PinitY   = Pinit.Y();
    PinitZ   = Pinit.Z();
    size=particle.size();

    
    os.str("");
    os_remove_g.str("");
    string pname[size];
    for(int i=0;i<size;i++){
      NucDeExParticle p = particle.at(i);
      PDG[i]=p._PDG;
      mass[i]=p._mass;
      totalE[i]=p.totalE();
      kE[i]=p.kE();
      PMag[i]=p._momentum.Mag();
      PX[i]=p._momentum.X();
      PY[i]=p._momentum.Y();
      PZ[i]=p._momentum.Z();
      flag[i]=p._flag;
      Ex_daughter[i]=p._Ex;
      // for decay mode string
      if(!p._flag) continue;  // intermediate state 
      if(p._name.length()>4) pname[i] = p._name.substr(0,1);
      else pname[i] = p._name;
      os << pname[i].c_str();
      if(i!=0 && pname[i] == "g") continue; // remove g if it is second
      os_remove_g << pname[i].c_str();
    }
    decay = os.str();
    decay_remove_g = os_remove_g.str();
    tree->Fill();

    // prepare detail fig
    if(flag_fig){
      TArrow* lXY[size], *lYZ[size], *lXZ[size];
      int color[size]={0};
      for(int i=0;i<size;i++){
        lXY[i] = new TArrow(0,0,PX[i],PY[i],0.01,"|>");
        lYZ[i] = new TArrow(0,0,PY[i],PZ[i],0.01,"|>");
        lXZ[i] = new TArrow(0,0,PX[i],PZ[i],0.01,"|>");
        lXY[i]->SetLineWidth(1);
        lYZ[i]->SetLineWidth(1);
        lXZ[i]->SetLineWidth(1);
        color[i]=400+1;
        for(int p=0;p<NucDeEx::num_particle;p++){
          if(PDG[i] == NucDeEx::PDG_particle[p]) color[i]=NucDeEx::color_root[p];
        }
        lXY[i]->SetLineColor(color[i]);
        lYZ[i]->SetLineColor(color[i]);
        lXZ[i]->SetLineColor(color[i]);
        lXY[i]->SetFillColor(color[i]);
        lYZ[i]->SetFillColor(color[i]);
        lXZ[i]->SetFillColor(color[i]);
      }
      TArrow* linitXY = new TArrow(0,0,PinitX,PinitY,0.015,"|>");
      TArrow* linitYZ = new TArrow(0,0,PinitY,PinitZ,0.015,"|>");
      TArrow* linitXZ = new TArrow(0,0,PinitX,PinitZ,0.015,"|>");
      linitXY->SetLineWidth(2);
      linitYZ->SetLineWidth(2);
      linitXZ->SetLineWidth(2);
      c_detail->Divide(2,2);
      //
      c_detail->cd(1);
      TH1F* waku1= gPad->DrawFrame(-200,-200,200,200);
      waku1->SetTitle("XY plane");
      waku1->GetXaxis()->SetTitle("Px (MeV)");
      waku1->GetYaxis()->SetTitle("Py (MeV)");
      linitXY->Draw();
      for(int i=0;i<size;i++){
        if(flag[i]) lXY[i]->Draw();
      }
      //
      c_detail->cd(2);
      TH1F* waku2= gPad->DrawFrame(-200,-200,200,200);
      waku2->SetTitle("YZ plane");
      waku2->GetXaxis()->SetTitle("Py (MeV)");
      waku2->GetYaxis()->SetTitle("Pz (MeV)");
      linitYZ->Draw();
      for(int i=0;i<size;i++){
        if(flag[i]) lYZ[i]->Draw();
      }
      //
      c_detail->cd(3);
      TH1F* waku3= gPad->DrawFrame(-200,-200,200,200);
      waku3->SetTitle("XZ plane");
      waku3->GetXaxis()->SetTitle("Px (MeV)");
      waku3->GetYaxis()->SetTitle("Pz (MeV)");
      linitXZ->Draw();
      for(int i=0;i<size;i++){
        if(flag[i]) lXZ[i]->Draw();
      }
      //
      c_detail->cd(4);
      TPaveText* t = new TPaveText(0.1,0.1,0.9,0.9);
      os.str("");
      os << "eventID = " << eventID;
      t->AddText(os.str().c_str());
      os.str("");
      os << argv[1] << ",  Ex = " << fixed << setprecision(2) << Ex
         << ",  shell = " << shell;
      t->AddText(os.str().c_str());
      os.str("");
      os << "P=(" << setprecision(1) << PinitX
         << ", " << setprecision(1) << PinitY
         << ", " << setprecision(1) << PinitZ << ")";
      t->AddText(os.str().c_str());
      //
      t->AddText("--- Decay products ---");
      for(int i=0;i<size;i++){
        if(!flag[i]) continue;
        os.str("");
        os << pname[i].c_str();
        os << ": kE = " <<  kE[i];
        //t->AddText(os.str().c_str());
        //((TText*)t->GetListOfLines()->Last())->SetTextColor(color[i]);
        os << ", P=(" << setprecision(1) << PX[i]
           << ", " << setprecision(1) << PY[i]
           << ", " << setprecision(1) << PZ[i] << ")";
          t->AddText(os.str().c_str());
        ((TText*)t->GetListOfLines()->Last())->SetTextColor(color[i]);
      }
      t->AddText("---------------------");
      t->AddText("LAB frame, MeV");
      ((TText*)t->GetListOfLines()->Last())->SetTextSize(0.045);
      ((TText*)t->GetListOfLines()->Last())->SetTextColor(kGray+1);
      t->Draw("same");
      //
      c_detail->Print(pdfname.c_str());
      c_detail->Update();
      c_detail->Clear();
    }
    eventID++;
  }
  if(flag_fig){
    c_detail->Print( (pdfname + (string)"]").c_str() );
    delete c_detail;
  }




  outf->cd();
  tree->Write();
  outf->Close();
  delete outf;



  return 0;
}