#include <iostream>
#include <string>

#include"event1.h"

#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace std;

int main(int argc, char* argv[]){
	// --------------------- //
	string prefix = "SF_LFG_14_1000_CCQE_C";
	if(argc==2) prefix = argv[1];
	// --------------------- //
  ostringstream os;
  os.str("");
	os << "output/" << prefix.c_str() << ".root";

  TFile* rootf = new TFile(os.str().c_str(),"READ");
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

	TH1D* h_Ev = new TH1D("h_Ev","",500,0,1e4);
	TH1D* h_Pinit = new TH1D("h_Pinit","",100,0,500);
  TH1D* h_PID_postFSI = new TH1D("h_PID_postFSI","",5000,-2500,2500);
  //TH1D* h_nmulti_preFSI = new TH1D("h_nmulti_preFSI","",10,-0.5,9.5);
  TH1D* h_nmulti_postFSI = new TH1D("h_nmulti_postFSI","",10,-0.5,9.5);
  TH1D* h_kEn = new TH1D("h_kEn","",500,0,1000);
  TH1D* h_kEp = new TH1D("h_kEp","",500,0,1000);
  TH1D* h_kEg = new TH1D("h_kEg","",200,0,10);
  TH1D* h_lowkEn = new TH1D("h_lowkEn","",100,0,50);
  TH1D* h_lowkEp = new TH1D("h_lowkEp","",100,0,50);
  TH1D* h_MissE = new TH1D("h_MissE","",200,0,200);
  TH1D* h_MissE_postFSI = new TH1D("h_MissE_postFSI","",200,0,200);
  TH1D* h_MissE_FSI = new TH1D("h_MissE_FSI","",400,-200,200);


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

		h_Ev->Fill(Ev);
		h_Pinit->Fill(event->in[1].momentum());

		h_nmulti_postFSI->Fill(event->fof(2112));

		// postFSI
		double Enuccsum=0;
		bool flag_pi=0;
		bool flag_src=0;
		for(int j=1;j<event->f();j++){ // w/o charged lepton
			int PDG = event->post[j].pdg;
			h_PID_postFSI->Fill(PDG);
			if(PDG==2112){
				h_kEn->Fill(event->post[j].Ek());
				h_lowkEn->Fill(event->post[j].Ek());
			}else if(PDG==2212){
				h_kEp->Fill(event->post[j].Ek());
				h_lowkEp->Fill(event->post[j].Ek());
			}
			Enuccsum += event->post[j].Ek();
			cout << "postFSI " << PDG << "   " << event->post[j].Ek()
					 << "    " << event->post[j].momentum()
					 << "   (" << event->post[j].p().x 
					 << "   " << event->post[j].p().y
					 << "   " << event->post[j].p().z << ")" << endl;
			if(abs(PDG)==211 || abs(PDG)==111) flag_pi=1;
			if(event->post[j].momentum()==event->in[1].momentum()) flag_src=1;
    }
		if(flag_pi || flag_src) continue;
		double MissE_postFSI = Ev-Evv-Enuccsum;
		h_MissE_postFSI->Fill(MissE_postFSI);

		cout << "el  " << event->number_of_nucleon_elastic() << endl;
		cout << "ce  " << event->number_of_nucleon_ce() << endl;
		cout << "spp  " << event->number_of_nucleon_spp() << endl;
		cout << "dpp  " << event->number_of_nucleon_dpp() << endl;
		
		if(event->f()>2){
			double MissE_FSI = event->out[1].Ek() - Enuccsum;		
			if(event->number_of_nucleon_elastic() + event->number_of_nucleon_ce() 
					+ event->number_of_nucleon_spp() + event->number_of_nucleon_dpp()>0){
				h_MissE_FSI->Fill(MissE_FSI);
				cout << MissE_FSI << endl;
			}
		}
		//cout << "W   " << event->W() << endl;

		cout << endl;
  }





	TCanvas* c_Ev = new TCanvas("c_Ev","c_Ev",0,0,800,600);
	h_Ev->GetXaxis()->SetTitle("Incoming neutrino energy (MeV)");
	h_Ev->GetYaxis()->SetTitle("Events/bin");
	h_Ev->GetYaxis()->SetMaxDigits(2);
	h_Ev->Draw("HIST");
	os.str("");
	os << "fig/fig_Ev_" << prefix.c_str() << ".pdf";
	c_Ev->Print(os.str().c_str());

	TCanvas* c_Pinit = new TCanvas("c_Pinit","c_Pinit",0,0,800,600);
	h_Pinit->GetXaxis()->SetTitle("Momentum of target nucleon (MeV)");
	h_Pinit->GetYaxis()->SetTitle("Events/bin");
	h_Pinit->GetYaxis()->SetMaxDigits(2);
	h_Pinit->Draw("HIST");
	os.str("");
	os << "fig/fig_Pinit_" << prefix.c_str() << ".pdf";
	c_Pinit->Print(os.str().c_str());


	TCanvas* c_PID_postFSI = new TCanvas("c_PID_postFSI","c_PID_postFSI",0,0,800,600);
	gPad->SetLogy();
	h_PID_postFSI->GetXaxis()->SetTitle("PID code (postFSI)");
	h_PID_postFSI->GetYaxis()->SetTitle("Particles/bin");
	h_PID_postFSI->Draw("HIST");
	os.str("");
	os << "fig/fig_PID_postFSI_" << prefix.c_str() << ".pdf";
	c_PID_postFSI->Print(os.str().c_str());


	TCanvas* c_nmulti_postFSI = new TCanvas("c_nmulti_postFSI","c_nmulti_postFSI",0,0,800,600);
	h_nmulti_postFSI->GetXaxis()->SetTitle("Neutron multiplicity");
	h_nmulti_postFSI->GetYaxis()->SetTitle("Events/bin");
	h_nmulti_postFSI->GetYaxis()->SetMaxDigits(2);
	h_nmulti_postFSI->Draw("HIST");
	os.str("");
	os << "Mean n multi. = " << fixed << setprecision(3) << h_nmulti_postFSI->GetMean();
	TLatex* l_mean_nmulti = new TLatex(3,h_nmulti_postFSI->GetMaximum()*0.7,os.str().c_str());
	l_mean_nmulti->Draw("same");
	os.str("");
	os << "fig/fig_nmulti_postFSI_" << prefix.c_str() << ".pdf";
	c_nmulti_postFSI->Print(os.str().c_str());



	TCanvas* c_kE = new TCanvas("c_kE","c_kE",0,0,1600,600);
	c_kE->Divide(2);
	c_kE->cd(1);
	h_kEp->SetTitle("Protons");
	h_kEp->GetXaxis()->SetTitle("Kinetic energy (MeV)");
	h_kEp->GetYaxis()->SetTitle("Particles/bin");
	h_kEp->GetYaxis()->SetMaxDigits(2);
	h_kEp->GetXaxis()->SetRangeUser(0,500);
	h_kEp->Draw("HIST");
	c_kE->cd(2);
	h_kEn->SetTitle("Neutrons");
	h_kEn->GetXaxis()->SetTitle("Kinetic energy (MeV)");
	h_kEn->GetYaxis()->SetTitle("Particles/bin");
	h_kEn->GetYaxis()->SetMaxDigits(2);
	h_kEn->GetXaxis()->SetRangeUser(0,500);
	h_kEn->Draw("HIST");
	os.str("");
	os << "fig/fig_kE_" << prefix.c_str() << ".pdf";
	c_kE->Print(os.str().c_str());

	if(h_kEg->GetEntries()>0){
		TCanvas* c_kEg = new TCanvas("c_kEg","c_kEg",0,0,800,600);
		h_kEg->SetTitle("Gammas");
		h_kEg->GetXaxis()->SetTitle("Kinetic energy (MeV)");
		h_kEg->GetYaxis()->SetTitle("Particles/bin");
		h_kEg->GetYaxis()->SetMaxDigits(2);
		h_kEg->Draw("HIST");
		os.str("");
		os << "fig/fig_kEg_" << prefix.c_str() << ".pdf";
		c_kEg->Print(os.str().c_str());
	}



	TCanvas* c_lowkE = new TCanvas("c_lowkE","c_lowkE",0,0,1600,600);
	c_lowkE->Divide(2);
	c_lowkE->cd(1);
	h_lowkEp->SetTitle("Protons");
	h_lowkEp->GetXaxis()->SetTitle("Kinetic energy (MeV)");
	h_lowkEp->GetYaxis()->SetTitle("Particles/bin");
	h_lowkEp->GetYaxis()->SetMaxDigits(2);
	h_lowkEp->Draw("HIST");
	c_lowkE->cd(2);
	h_lowkEn->SetTitle("Neutrons");
	h_lowkEn->GetXaxis()->SetTitle("Kinetic energy (MeV)");
	h_lowkEn->GetYaxis()->SetTitle("Particles/bin");
	h_lowkEn->GetYaxis()->SetMaxDigits(2);
	h_lowkEn->Draw("HIST");
	os.str("");
	os << "fig/fig_lowkE_" << prefix.c_str() << ".pdf";
	c_lowkE->Print(os.str().c_str());



	TCanvas* c_MissE = new TCanvas("c_MissE","c_MissE",0,0,800,600);
	h_MissE->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_MissE->GetYaxis()->SetTitle("Events/bin");
	h_MissE->GetYaxis()->SetMaxDigits(2);
	h_MissE->Draw("HIST");
	os.str("");
	os << "fig/fig_MissE_" << prefix.c_str() << ".pdf";
	c_MissE->Print(os.str().c_str());


	TCanvas* c_MissE_postFSI = new TCanvas("c_MissE_postFSI","c_MissE_postFSI",0,0,800,600);
	h_MissE_postFSI->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_MissE_postFSI->GetYaxis()->SetTitle("Events/bin");
	h_MissE_postFSI->GetYaxis()->SetMaxDigits(2);
	h_MissE_postFSI->Draw("HIST");
	os.str("");
	os << "fig/fig_MissE_postFSI_" << prefix.c_str() << ".pdf";
	c_MissE_postFSI->Print(os.str().c_str());

	TCanvas* c_MissE_FSI = new TCanvas("c_MissE_FSI","c_MissE_FSI",0,0,800,600);
	h_MissE_FSI->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_MissE_FSI->GetYaxis()->SetTitle("Events/bin");
	h_MissE_FSI->GetYaxis()->SetMaxDigits(2);
	h_MissE_FSI->Draw("HIST");
	os.str("");
	os << "fig/fig_MissE_FSI_" << prefix.c_str() << ".pdf";
	c_MissE_FSI->Print(os.str().c_str());



	os.str("");
	os << "output/histogram_" << prefix.c_str() << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	h_Pinit->Write();
	h_nmulti_postFSI->Write();
	h_kEn->Write();
	h_kEp->Write();
	h_lowkEn->Write();
	h_lowkEp->Write();
	h_MissE->Write();
	h_MissE_postFSI->Write();
	h_MissE_FSI->Write();
	outf->Close();
	delete outf;

	rootf->Close();
	delete rootf;

  return 0;
}
