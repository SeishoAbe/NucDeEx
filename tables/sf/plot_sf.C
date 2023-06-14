#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int plot_sf(){
  int flag=0;
  // 0: 12C
  // 1: 16O
  
  string nucleus;
  int A;
  if(flag==0){
    nucleus = "C";
    A=12;
  }else if(flag==1){
    nucleus = "O";
    A=16;
  }else{
    return 1;
  }


  ostringstream os;
  os.str("");
  os << "pke" << A;
  if(flag==0) os << "_tot";
	string prefix = os.str();
	os << ".grid";
  ifstream ifs(os.str().c_str());
  if(!ifs.is_open()){
    cerr << "Cannot open " << os.str().c_str() << endl;
    return 1;
  }
  char buf[256];

  
  int bin_E, bin_p;
  double min_E, min_p;
  double max_E, max_p;
  string st;
  // 
  ifs.getline(buf,sizeof(buf));
  st = (string) buf;
  istringstream(st) >> bin_E >> bin_p;
  // 
  ifs.getline(buf,sizeof(buf));
  st = (string) buf;
  istringstream(st) >> min_E >> min_p;
  // 
  ifs.getline(buf,sizeof(buf));
  st = (string) buf;
  istringstream(st) >> max_E >> max_p;


  cout << "Energy: (" << bin_E << ", " << min_E << ", " << max_E << ")" << endl;
  cout << "Momentum: (" << bin_p << ", " << min_p << ", " << max_p << ")" << endl;

  TH2D* h_sf_nominal = new TH2D("h_sf_nominal","",bin_p,min_p,max_p,bin_E,min_E,max_E);
      // vector p
  TH2D* h_sf = new TH2D("h_sf","",bin_p,min_p,max_p,bin_E,min_E,max_E); 
      // scalar p
  
  TH1D* h_sf_p = new TH1D("h_sf_p","",bin_p,min_p,max_p);
  TH1D* h_sf_E = new TH1D("h_sf_E","",bin_E,min_E,max_E);
  
  int line = bin_E/4;
  if(flag==1) line++;
  for(int i=0;i<bin_p;i++){
    double p;
    ifs.getline(buf,sizeof(buf));
    st = (string) buf;
    istringstream(st) >> p; // momentum
    cout << p << endl;
    for(int i=0;i<line;i++){
      ifs.getline(buf,sizeof(buf));
      double e1,e2,e3,e4;
      double prob1,prob2,prob3,prob4;
      st = (string) buf;
      if(i==bin_E/4 && flag==1){
        istringstream(st) >> e1 >> prob1
                          >> e2 >> prob2;
        h_sf_nominal->Fill(p,e1,prob1);
        h_sf_nominal->Fill(p,e2,prob2);
        //
        h_sf->Fill(p,e1,prob1*pow(p,2));
        h_sf->Fill(p,e2,prob2*pow(p,2));
        //
        h_sf_p->Fill(p,prob1*pow(p,2));
        h_sf_p->Fill(p,prob2*pow(p,2));
        //
        h_sf_E->Fill(e1,prob1*pow(p,2));
        h_sf_E->Fill(e2,prob2*pow(p,2));
      }else{
        istringstream(st) >> e1 >> prob1
                          >> e2 >> prob2
                          >> e3 >> prob3
                          >> e4 >> prob4;
        h_sf_nominal->Fill(p,e1,prob1);
        h_sf_nominal->Fill(p,e2,prob2);
        h_sf_nominal->Fill(p,e3,prob3);
        h_sf_nominal->Fill(p,e4,prob4);
        //
        h_sf->Fill(p,e1,prob1*pow(p,2));
        h_sf->Fill(p,e2,prob2*pow(p,2));
        h_sf->Fill(p,e3,prob3*pow(p,2));
        h_sf->Fill(p,e4,prob4*pow(p,2));
        //
        h_sf_p->Fill(p,prob1*pow(p,2));
        h_sf_p->Fill(p,prob2*pow(p,2));
        h_sf_p->Fill(p,prob3*pow(p,2));
        h_sf_p->Fill(p,prob4*pow(p,2));
        //
        h_sf_E->Fill(e1,prob1*pow(p,2));
        h_sf_E->Fill(e2,prob2*pow(p,2));
        h_sf_E->Fill(e3,prob3*pow(p,2));
        h_sf_E->Fill(e4,prob4*pow(p,2));
      }
    }
  }

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  TCanvas* c_nominal = new TCanvas("c_nominal","c_nominal",0,0,800,600);
  gPad->SetLogz();
  h_sf_nominal->SetStats(0);
  h_sf_nominal->GetXaxis()->SetTitle("Momentum (vector)(MeV)");
  h_sf_nominal->GetYaxis()->SetTitle("Removal energy (MeV)");
  h_sf_nominal->Draw("colz");
  os.str("");
  os << "fig_sf_nominal_" << nucleus.c_str() << ".pdf";
  c_nominal->Update();
  c_nominal->Print(os.str().c_str());

  TCanvas* c = new TCanvas("c","c",0,0,800,600);
  gPad->SetLogz();
  h_sf->SetStats(0);
  h_sf->GetXaxis()->SetTitle("Momentum (MeV)");
  h_sf->GetYaxis()->SetTitle("Removal energy (MeV)");
  h_sf->Draw("colz");
  os.str("");
  os << "fig_sf_" << nucleus.c_str() << ".pdf";
  c->Update();
  c->Print(os.str().c_str());


  TCanvas* c_1D = new TCanvas("c_1D","c_1D",0,0,1200,600);
  c_1D->Divide(2);
  c_1D->cd(1);
  h_sf_p->GetXaxis()->SetTitle("Momentum (MeV)");
  h_sf_p->GetYaxis()->SetTitle("A.U.");
  h_sf_p->SetStats(0);
  h_sf_p->Draw("HIST");
  c_1D->cd(2);
  h_sf_E->GetXaxis()->SetTitle("Removal energy (MeV)");
  h_sf_E->GetYaxis()->SetTitle("A.U.");
  h_sf_E->SetStats(0);
  h_sf_E->GetXaxis()->SetRangeUser(0,100);
  h_sf_E->Draw("HIST");
  os.str("");
  os << "fig_sf_1D_" << nucleus.c_str() << ".pdf";
  c_1D->Update();
  c_1D->Print(os.str().c_str());

  os.str("");
  os << prefix.c_str() << ".root";
  TFile* outf = new TFile(os.str().c_str(),"RECREATE");
  h_sf_nominal->Write();
  h_sf->Write();
  h_sf_p->Write();
  h_sf_E->Write();
  outf->Close();
  delete outf;


  return 0;
}
