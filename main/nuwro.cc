#include <iostream>
#include <string>

#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"

#include"event1.h"

#include "NucleusTable.hh"
#include "Deexcitation.hh"

using namespace std;

int main(int argc, char* argv[]){
	//------------------//
	string prefix = "SF_LFG_14_1000_CCQE_C";
	if(argc==2) prefix = argv[1];
	int seed=1;
	//-----------------.//

	ostringstream os;

	// --- prepare deex tools
	int Zt, Nt;
	double S;
	Deexcitation* deex = new Deexcitation(2, 1);
	deex->SetSeed(seed);
	deex->SetVerbose(1);
  NucleusTable* nucleus_table = deex->GetNucleusTablePtr();
	if(prefix.find("_C_")!=string::npos){
		Zt=6;
		Nt=6;
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke12_tot.root";
		S = nucleus_table->GetNucleusPtr("12C")->S[2];
	}else if(prefix.find("_O_")!=string::npos){
		Zt=8;
		Nt=8;
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke16.root";
		S = nucleus_table->GetNucleusPtr("16O")->S[2];
	}
	cout << "Separation E = " << S << endl;
	TFile* root = new TFile(os.str().c_str(),"READ");
	TH2D* h_sf_int = (TH2D*) root->Get("h_sf_int");
	h_sf_int->SetDirectory(0);
	gRandom->SetSeed(seed); // for TH2 GetRandom2
	root->Close();
	delete root;

	// --- input neut root
	os.str("");
	os << "output_nuwro/" << prefix.c_str() << ".root";
  TFile* rootf = new TFile(os.str().c_str(),"READ");
	cout << os.str().c_str() << endl;
  TTree* tree = (TTree*) rootf->Get("treeout");
  event* event=0;
  tree->SetBranchAddress("e",&event);

  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.08);
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetLabelSize(0.05,"XYZ");
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetLegendFont(132);
  gStyle->SetTitleYOffset(0.95);

	TH1D* h_nmulti_postFSI = new TH1D("h_nmulti_postFSI","",10,-0.5,9.5);
	TH1D* h_nmulti_postdeex = new TH1D("h_nmulti_postdeex","",10,-0.5,9.5);
	TH1D* h_MissE = new TH1D("h_MissE","",200,0,200);
	TH1D* h_Ex[4]; // single [0]: all
	for(int i=0;i<4;i++){
		os.str("");
		os << "h_Ex_" << i;
		h_Ex[i] = new TH1D(os.str().c_str(),"",500,-100,400);
	}
	TH1D* h_Ex_multi = new TH1D("h_Ex_multi","",500,-100,400); // multi nucleon hole

	// --- read root
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    // -- calculate missing E (pre-FSI)
    double Ev = event->in[0].E();
    double Evv = event->out[0].E(); // mass E inc
    double Enucc = event->out[1].E(); // mass E inc
		double massnuc=event->in[1].mass();
    double MissE = Ev-Evv-Enucc + massnuc;
		h_MissE->Fill(MissE);
		
		cout << "Init    " << event->in[1].pdg << "   " << event->in[1].Ek()
				 << "    " << event->in[1].momentum() 
				 << "   (" << event->in[1].p().x 
				 << "   " << event->in[1].p().y
				 << "   " << event->in[1].p().z << ")" << endl;
		cout << "preFSI  " << event->out[1].pdg << "   " << event->out[1].Ek()
		     << "    " << event->out[1].momentum()
				 << "   (" << event->out[1].p().x 
				 << "   " << event->out[1].p().y
				 << "   " << event->out[1].p().z << ")" << endl;

		h_nmulti_postFSI->Fill(event->fof(2112));
  }

	TCanvas* c_nmulti_postFSI = new TCanvas("c_nmulti_postFSI","c_nmulti_postFSI",0,0,800,600);
	h_nmulti_postFSI->GetXaxis()->SetTitle("Neutron multiplicity");
	h_nmulti_postFSI->GetYaxis()->SetTitle("Events/bin");
	h_nmulti_postFSI->GetYaxis()->SetMaxDigits(3);
	h_nmulti_postFSI->SetStats(0);
	h_nmulti_postFSI->Scale(1./h_nmulti_postFSI->GetEntries());
	h_nmulti_postFSI->Draw("HIST");
	h_nmulti_postdeex->SetLineColor(kRed);
	h_nmulti_postdeex->Scale(1./h_nmulti_postdeex->GetEntries());
	h_nmulti_postdeex->Draw("HISTsame");
	os.str("");
	os << "Mean n multi. = " << fixed << setprecision(3) << h_nmulti_postFSI->GetMean();
	TLatex* l_mean_nmulti = new TLatex(3,h_nmulti_postFSI->GetMaximum()*0.7,os.str().c_str());
	l_mean_nmulti->Draw("same");
	os.str("");
	os << "Mean n multi. = " << fixed << setprecision(3) << h_nmulti_postdeex->GetMean();
	TLatex* l_mean_nmulti_postdeex = new TLatex(3,h_nmulti_postFSI->GetMaximum()*0.6,os.str().c_str());
	l_mean_nmulti_postdeex->SetTextColor(kRed);
	l_mean_nmulti_postdeex->Draw("same");
	os.str("");
	os << "fig_nuwro/fig_nmulti_postFSI_" << prefix.c_str() << ".pdf";
	c_nmulti_postFSI->Print(os.str().c_str());

	TCanvas* c_MissE = new TCanvas("c_MissE","c_MissE",0,0,800,600);
	h_MissE->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_MissE->GetYaxis()->SetTitle("Events/bin");
	h_MissE->GetYaxis()->SetMaxDigits(2);
	h_MissE->Draw("HIST");
	os.str("");
	os << "fig_nuwro/fig_MissE_" << prefix.c_str() << ".pdf";
	c_MissE->Print(os.str().c_str());

	TCanvas* c_Ex = new TCanvas("c_Ex","c_Ex",0,0,800,600);
	h_Ex[0]->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_Ex[0]->GetYaxis()->SetTitle("Events/bin");
	h_Ex[0]->GetYaxis()->SetMaxDigits(2);
	h_Ex[0]->GetXaxis()->SetRangeUser(-10,100);
	h_Ex[0]->SetStats(0);
	h_Ex[0]->Draw("HIST");
	THStack* h_s_Ex = new THStack("h_s_Ex","");
	const int color[4]={1,600-7,632-7,920};
	for(int i=1;i<4;i++){
		h_Ex[i]->SetFillColor(color[i]);
		h_s_Ex->Add(h_Ex[i]);
	}
	h_s_Ex->Draw("same");
	//
	os.str("");
	os << "Prob(s1/2)=" << fixed << setprecision(1) << (double)h_Ex[1]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
	TText* t_s12 = new TText(40,h_Ex[0]->GetMaximum()*0.9,os.str().c_str());
	t_s12->Draw("same");
	//
	os.str("");
	os << "Prob(p3/2)=" << fixed << setprecision(1) << (double)h_Ex[2]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
	TText* t_p32 = new TText(40,h_Ex[0]->GetMaximum()*0.8,os.str().c_str());
	t_p32->Draw("same");
	if(prefix.find("_O_")!=string::npos){
		os.str("");
		os << "Prob(p1/2)=" << fixed << setprecision(1) << (double)h_Ex[3]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
		TText* t_p12 = new TText(40,h_Ex[0]->GetMaximum()*0.7,os.str().c_str());
		t_p12->Draw("same");
	}
	//
	gPad->RedrawAxis();
	os.str("");
	os << "fig_nuwro/fig_Ex_" << prefix.c_str() << ".pdf";
	c_Ex->Print(os.str().c_str());
	
	// --- save 
	os.str("");
	os << "output_nuwro/histogram_deex_" << prefix.c_str() << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	h_nmulti_postFSI->Write();
	h_nmulti_postdeex->Write();
	h_MissE->Write();
	for(int i=0;i<4;i++){
		h_Ex[i]->Write();
	}
	h_Ex_multi->Write();
	outf->Close();
	delete outf;

	rootf->Close();
	delete rootf;

	return 0;
}
