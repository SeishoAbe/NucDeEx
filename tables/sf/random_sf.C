#include "../../include/NucleusTable.hh"
#include "../../include/Deexcitation.hh"

R__LOAD_LIBRARY(../../lib/libTALYStool);

using namespace std;

int random_sf(){
	// --------------------//
	const int numofevent=1e6; // to be generated

	const double Ex_p32 =4.0;
	const double Ex_s12 =14;
	string target="16O";
	//
	/*

	const double Ex_p32 =-1e9;
	const double Ex_s12 =16.0;
	string target="12C";
	*/
	double max_mom=500;
	// --------------------//

	const int ldmodel=2;
	const bool parity_optmodall=1;
	int seed=1; // default: 1
	ostringstream os,os_remove_g;

	// Set deex tool
	Deexcitation* deex = new Deexcitation(ldmodel, parity_optmodall);
	deex->SetSeed(seed); // 0: time
	deex->SetVerbose(1);
  NucleusTable* nucleus_table = deex->GetNucleusTablePtr();
	Nucleus* nuc = nucleus_table->GetNucleusPtr(target.c_str());
	gStyle->SetTextSize(0.08);
	gStyle->SetTitleSize(0.045);
	gStyle->SetTitleXSize(0.045);
	gStyle->SetTitleYSize(0.045);
	gStyle->SetTitleYOffset(0.95);
	TH2D* h_sf_random = new TH2D("h_sf_random","",400,0,800,400,0,400);
	TH1D* h_sf_E_random = new TH1D("h_sf_E_random","",400,0,400); // missing E
	TH1D* h_sf_p_random = new TH1D("h_sf_p_random","",400,0,800);
	TH1D* h_sf_Ex_random = new TH1D("h_sf_Ex_random","",500,-100,400);  // excitation E

	// Set Ex and mom tables (from SF) and proton separation energy
	double MissE, Ex, S;
	os.str("");	
	if(target=="12C"){
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke12_tot.root";
		S = nucleus_table->GetNucleusPtr("12C")->S[2];
	}else if(target=="16O"){
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke16.root";
		S = nucleus_table->GetNucleusPtr("16O")->S[2];
	}
	cout << "S = " << S << endl;
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	TH2D* h_sf_int = (TH2D*) rootf->Get("h_sf_int");
	h_sf_int->SetDirectory(0);
	gRandom->SetSeed(seed); // for TH2 GetRandom2
	rootf->Close();
	delete rootf;


	//
	int eventID=0;
	TVector3 Pinit;
	double PinitMag;

	while(eventID<numofevent){
		// determine momentum (scalar) and missing E according to SF
		h_sf_int->GetRandom2(PinitMag,MissE);
		Ex=MissE-S;
		h_sf_random->Fill(PinitMag,MissE);
		h_sf_E_random->Fill(MissE);
		h_sf_p_random->Fill(PinitMag);
		h_sf_Ex_random->Fill(Ex);
		eventID++;
	}
	
	int total=0;
	int p12=0, p32=0, s12=0;
	for(int b=1;b<=500;b++){
		int cont = h_sf_Ex_random->GetBinContent(b);
		double Ex = h_sf_Ex_random->GetBinCenter(b);
		total += cont;
		if(Ex>Ex_s12) s12+=cont;
		else if(Ex>Ex_p32) p32+=cont;
		else p12+=cont;
	}
	cout << "p12 = " << p12 << endl;
	cout << "p32 = " << p32 << endl;
	cout << "s12 = " << s12 << endl;
	cout << "Prob(p12) = " << (double)p12/total*100 << endl;
	cout << "Prob(p32) = " << (double)p32/total*100 << endl;
	cout << "Prob(s12) = " << (double)s12/total*100 << endl;


	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	c->Divide(2,2);
	c->cd(1);
	h_sf_Ex_random->GetXaxis()->SetRangeUser(-10,100);
	h_sf_Ex_random->SetStats(0);
	h_sf_Ex_random->SetMinimum(0);
	h_sf_Ex_random->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_sf_Ex_random->GetYaxis()->SetTitle("Events/bin");
	h_sf_Ex_random->GetYaxis()->SetMaxDigits(2);
	h_sf_Ex_random->Draw("HIST");
	TLine* l_Ex_p32 = new TLine(Ex_p32,0,Ex_p32,h_sf_Ex_random->GetMaximum()*1.05);
	l_Ex_p32->SetLineStyle(2);
	l_Ex_p32->SetLineWidth(2);
	l_Ex_p32->SetLineColor(kRed);
	if(Ex_p32>0) l_Ex_p32->Draw("same");
	TLine* l_Ex_s12 = new TLine(Ex_s12,0,Ex_s12,h_sf_Ex_random->GetMaximum()*1.05);
	l_Ex_s12->SetLineStyle(2);
	l_Ex_s12->SetLineWidth(2);
	l_Ex_s12->SetLineColor(kRed);
	if(Ex_s12>0) l_Ex_s12->Draw("same");
	os.str("");
	os << "Prob(p3/2) = " << fixed << setprecision(1) << (double)p12/total*100 << "%";
	TText* t_p12 = new TText(40,h_sf_Ex_random->GetMaximum()*0.5,os.str().c_str());
	if(p12>0) t_p12->Draw("same");
	os.str("");
	os << "Prob(p3/2) = " << fixed << setprecision(1) << (double)p32/total*100 << "%";
	TText* t_p32 = new TText(40,h_sf_Ex_random->GetMaximum()*0.35,os.str().c_str());
	t_p32->Draw("same");
	os.str("");
	os << "Prob(s1/2) = " << fixed << setprecision(1) << (double)s12/total*100 << "%";
	TText* t_s12 = new TText(40,h_sf_Ex_random->GetMaximum()*0.2,os.str().c_str());
	t_s12->Draw("same");
	//
	os.str("");
	os << "Ex (p1/2<->p3/2) = " << Ex_p32 << " MeV";
	TText* t_Ex_p32 = new TText(20,h_sf_Ex_random->GetMaximum()*0.9,os.str().c_str());
	t_Ex_p32->SetTextColor(kRed);
	if(Ex_p32>0) t_Ex_p32->Draw("same");
	os.str("");
	os << "Ex (p3/2<->s1/2) = " << Ex_s12 << " MeV";
	TText* t_Ex_s12 = new TText(20,h_sf_Ex_random->GetMaximum()*0.7,os.str().c_str());
	t_Ex_s12->SetTextColor(kRed);
	if(Ex_s12>0) t_Ex_s12->Draw("same");
	//
	c->cd(2);
	gPad->SetLogz();
	h_sf_random->SetStats(0);
	h_sf_random->SetMinimum(1);
	h_sf_random->GetXaxis()->SetTitle("Momentum (MeV)");
	h_sf_random->GetYaxis()->SetTitle("Missing energy (MeV)");
	h_sf_random->GetXaxis()->SetRangeUser(0,max_mom);
	h_sf_random->GetYaxis()->SetRangeUser(0,100);
	h_sf_random->Draw("colz");
	h_sf_random->SetStats(0);
	TLine* l_sf_s12 = new TLine(0,Ex_s12+S,max_mom,Ex_s12+S);
	l_sf_s12->SetLineStyle(2);
	l_sf_s12->SetLineWidth(2);
	l_sf_s12->SetLineColor(kRed);
	if(Ex_s12>0) l_sf_s12->Draw("same");
	TLine* l_sf_p32 = new TLine(0,Ex_p32+S,max_mom,Ex_p32+S);
	l_sf_p32->SetLineStyle(2);
	l_sf_p32->SetLineWidth(2);
	l_sf_p32->SetLineColor(kRed);
	if(Ex_p32>0) l_sf_p32->Draw("same");
	c->cd(3);
	h_sf_E_random->SetStats(0);
	h_sf_E_random->SetMinimum(0);
	h_sf_E_random->GetXaxis()->SetRangeUser(0,100);
	h_sf_E_random->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_sf_E_random->GetYaxis()->SetTitle("Events/bin");
	h_sf_E_random->Draw("HIST");
	TLine* l_S = new TLine(S,0,S,h_sf_E_random->GetMaximum()*1.05);
	l_S->SetLineStyle(2);
	l_S->SetLineWidth(2);
	l_S->SetLineColor(kBlue);
	l_S->Draw("same");
	os.str("");
	os << "S(p) = " << S << " MeV";
	TText* t_S = new TText(40,h_sf_E_random->GetMaximum()*0.9,os.str().c_str());
	t_S->SetTextColor(kBlue);
	t_S->Draw("same");
	c->cd(4);
	h_sf_p_random->SetStats(0);
	h_sf_p_random->SetMinimum(0);
	h_sf_p_random->GetXaxis()->SetTitle("Momentum (MeV)");
	h_sf_p_random->GetYaxis()->SetTitle("Events/bin");
	h_sf_p_random->Draw("HIST");
	os.str("");
	os << "figure/fig_sf_random_" << target.c_str() << ".pdf";
	c->Print(os.str().c_str());

	return 0;
}
