#include <iostream>
#include <fstream>
#include <string>

using namespace std;

//double FUNC(TGraph* g,double*x, double*){
TGraph* g_sf_p, *g_sf_E;
double FUNC_p(double*x, double*)
{
	return g_sf_p->Eval(x[0]);
	//return g_sf_p->Eval(x[0],0,"S");
}
double FUNC_E(double*x, double*)
{
	return g_sf_E->Eval(x[0]);
	//return g_sf_E->Eval(x[0],0,"S");
}

TH2D* h_sf;
double FUNC_SF(double *val,double *par)
{
	double X = val[0];
	double Y = val[1];
	//int binX = h_sf->GetXaxis()->FindBin(X);
	//int binY = h_sf->GetYaxis()->FindBin(Y);
	//double Z = h_sf->GetBinContent(binX,binY);
	double Z = h_sf->Interpolate(X,Y);
	if(Y<15) Z=0;
	return Z;
}

TGraph2D* g_sf;
double FUNC_SF_G(double *val,double *par)
{
	double X = val[0];
	double Y = val[1];
	double Z = g_sf->Interpolate(X,Y);
	return Z;
	//return sqrt(X*X+Y*Y);
}

int plot_sf(){
  int flag=1;
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
  h_sf = new TH2D("h_sf","",bin_p,min_p,max_p,bin_E,min_E,max_E); 
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
    //cout << p << endl;
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

	g_sf_p = new TGraph();
	g_sf_E = new TGraph();
	g_sf_p->SetName("g_sf_p");
	g_sf_E->SetName("g_sf_E");
	g_sf_p->SetPoint(0,0,0);
	g_sf_E->SetPoint(0,0,0);
	for(int i=1;i<=bin_p;i++){
		g_sf_p->SetPoint(i,h_sf_p->GetBinCenter(i),h_sf_p->GetBinContent(i));
	}
	for(int i=1;i<=bin_E;i++){
		g_sf_E->SetPoint(i,h_sf_E->GetBinCenter(i),h_sf_E->GetBinContent(i));
	}

	g_sf = new TGraph2D();
	int index=0;
	g_sf->SetPoint(index,0,0,0);
	index++;
	double x,y;
	for(int j=1;j<=bin_E;j++){
		y = h_sf->GetYaxis()->GetBinCenter(j);
		g_sf->SetPoint(index,0,y,0);
		index++;
	}
	for(int i=1;i<=bin_p;i++){
		x = h_sf->GetXaxis()->GetBinCenter(i);
		g_sf->SetPoint(index,x,0,0);
		index++;
		for(int j=1;j<=bin_E;j++){
			y = h_sf->GetYaxis()->GetBinCenter(j);
			g_sf->SetPoint(index,x,y,h_sf->GetBinContent(i,j));
			index++;
		}
	}


	TF1* f_sf_p = new TF1("f_sf_p",FUNC_p,min_p,max_p,0);
	TF1* f_sf_E = new TF1("f_sf_E",FUNC_E,min_E,max_E,0);
	TF2* f_sf = new TF2("f_sf",FUNC_SF,min_p,max_p,min_E,max_E); // use th2
	//TF2* f_sf = new TF2("f_sf",FUNC_SF_G,min_p,max_p,min_E,max_E); // use tgraph2d
	
	int bin_p_int = bin_p*4;
	int bin_E_int = bin_E*5;
  cout << "INT_Energy: (" << bin_E_int << ", " << min_E << ", " << max_E << ")" << endl;
  cout << "INT_Momentum: (" << bin_p_int << ", " << min_p << ", " << max_p << ")" << endl;
  TH2D* h_sf_int = new TH2D("h_sf_int","",bin_p_int,min_p,max_p,bin_E_int,min_E,max_E); 
	for(int i=1;i<bin_p_int;i++){
		for(int j=1;j<bin_E_int;j++){
			double X,Y;
			X=h_sf_int->GetXaxis()->GetBinCenter(i);
			Y=h_sf_int->GetYaxis()->GetBinCenter(j);
			h_sf_int->SetBinContent(i,j,h_sf->Interpolate(X,Y));
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
  os << "figure/fig_sf_nominal_" << nucleus.c_str() << ".pdf";
  c_nominal->Update();
  c_nominal->Print(os.str().c_str());

  TCanvas* c = new TCanvas("c","c",0,0,800,600);
  gPad->SetLogz();
  h_sf->SetStats(0);
  h_sf->GetXaxis()->SetTitle("Momentum (MeV)");
  h_sf->GetYaxis()->SetTitle("Removal energy (MeV)");
  h_sf->Draw("colz");
	g_sf->SetMarkerStyle(8);
	g_sf->SetMarkerColor(kBlack);
	//g_sf->Draw("PCOL");
	//g_sf->Draw("surf");
	f_sf->Draw("same");
  os.str("");
  os << "figure/fig_sf_" << nucleus.c_str() << ".pdf";
  c->Update();
  c->Print(os.str().c_str());

	TCanvas* c_sf_g2D = new TCanvas("c_sf_g_2D","",0,0,800,600);
	gPad->SetLogz();
	//g_sf->GetYaxis()->SetRangeUser(0,50);
	g_sf->Draw("SURF1");

	
  TCanvas* c_1D = new TCanvas("c_1D","c_1D",0,0,1200,600);
  c_1D->Divide(2);
  c_1D->cd(1);
  h_sf_p->GetXaxis()->SetTitle("Momentum (MeV)");
  h_sf_p->GetYaxis()->SetTitle("A.U.");
  h_sf_p->SetStats(0);
  h_sf_p->Draw("HIST");
	g_sf_p->SetLineColor(kRed);
	g_sf_p->SetMarkerColor(kRed);
	g_sf_p->SetMarkerStyle(8);
	g_sf_p->Draw("PLsame");
	f_sf_p->SetLineColor(kBlue);
	f_sf_p->Draw("same");
  c_1D->cd(2);
  h_sf_E->GetXaxis()->SetTitle("Removal energy (MeV)");
  h_sf_E->GetYaxis()->SetTitle("A.U.");
  h_sf_E->SetStats(0);
  h_sf_E->GetXaxis()->SetRangeUser(0,100);
  h_sf_E->Draw("HIST");
	g_sf_E->SetLineColor(kRed);
	g_sf_E->SetMarkerColor(kRed);
	g_sf_E->SetMarkerStyle(8);
	g_sf_E->Draw("PLsame");
	f_sf_E->SetLineColor(kBlue);
	f_sf_E->Draw("same");
  os.str("");
  os << "figure/fig_sf_sf_1D_" << nucleus.c_str() << ".pdf";
  c_1D->Update();
  c_1D->Print(os.str().c_str());
	

	TCanvas* c_int = new TCanvas("c_int","c_int",0,0,800,600);
	gPad->SetLogz();
	h_sf_int->SetStats(0);
	h_sf_int->Draw("colz");
	


  os.str("");
  os << prefix.c_str() << ".root";
  TFile* outf = new TFile(os.str().c_str(),"RECREATE");
  h_sf_nominal->Write();
  h_sf->Write();
  h_sf_p->Write();
  h_sf_E->Write();
	g_sf_p->Write();
	g_sf_E->Write();
	f_sf_p->Write();
	f_sf_E->Write();
	f_sf->Write();
	h_sf_int->Write();
  outf->Close();
  //delete outf;




  return 0;
}
