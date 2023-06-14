#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "NucleusTable.hh"
#include "Deexcitation.hh"

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH2.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TLine.h>

using namespace std;

int main(int argc, char* argv[]){
	if(argc<=1){
		cerr << "Input: " << argv[0] << " [Target nucleus] [Random Seed (optional)]" << endl;
		return 0;
	}
	int seed=1; // default: 1
	if(argc==3) seed = atoi(argv[2]);

	// ---- FIXME --- // 
	const int numofevent=1e4; // to be generated
	const double Ex_min = 16.;
	const double Ex_max = 35.;
	// -------------- //

	ostringstream os;

	// Get Z and N
  NucleusTable* nucleus_table = new NucleusTable();
  if(!nucleus_table->ReadTables()){
		cerr << "something wrong" << endl;
		return 1;
	}
	Nucleus* nuc = nucleus_table->GetNucleusPtr(argv[1]);
	const int Z = nuc->Z;
	const int N = nuc->N;
	
	// Set deex tool
	Deexcitation* deex = new Deexcitation();
	deex->SetSeed(seed); // 0: time
	deex->SetVerbose(1);
	cout << "SEED = " << seed << endl;
	


	// prepare output root file
	os.str("");	
	os << "sim_out/output_" << argv[1] << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	TTree* tree = new TTree("tree",""); // LAB Freme, MeV
	int eventID, size;
	double MissE, Ex, S;
	double PinitMag, PinitX,PinitY,PinitZ;
	//
	int PDG[bins];
	double mass[bins];
	double totalE[bins],kE[bins];
	double PMag[bins], PX[bins],PY[bins],PZ[bins];
	tree->Branch("eventID",&eventID,"eventID/I");
	tree->Branch("MissE",&MissE,"MissE/D");
	tree->Branch("S",&S,"S/D");
	tree->Branch("Ex",&Ex,"Ex/D");
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

	TH2D* h_sf_random = new TH2D("h_sf_random","",800,0,800,400,0,400);
	TH1D* h_sf_E_random = new TH1D("h_sf_E_random","",400,0,400); // missing E
	TH1D* h_sf_p_random = new TH1D("h_sf_p_random","",800,0,800);
	TH1D* h_sf_Ex_random = new TH1D("h_sf_Ex_random","",400,0,400);  // excitation E


	// Set Ex and mom tables (sf)
	os.str();
	if(Z+N==11){ // 12C
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke12_tot.root";
		S = nucleus_table->GetNucleusPtr("12C")->S[2];
	}
	if(Z+N==15){
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke16.root";
		S = nucleus_table->GetNucleusPtr("16O")->S[2];
	}
	cout << "S = " << S << endl;
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	TH2D* h_sf_int = (TH2D*) rootf->Get("h_sf_int");
	cout << gRandom->GetSeed() << endl;
	gRandom->SetSeed(seed); // for TH2 GetRandom2


	
	while(deex->GetEventID()<numofevent){
		// determine momentum (scalar) and missing E according to SF
		h_sf_int->GetRandom2(PinitMag,MissE);
		Ex=MissE-S;
		h_sf_random->Fill(PinitMag,MissE);
		h_sf_E_random->Fill(MissE);
		h_sf_p_random->Fill(PinitMag);
		h_sf_Ex_random->Fill(Ex);
		// determine angle 
		double costheta = 2.*gRandom->Rndm()-1;
		double sintheta = sqrt( 1. - pow(costheta,2) );
		double phi      = 2*TMath::Pi()*gRandom->Rndm();
		TVector3 Pinit(PinitMag*sintheta*cos(phi),
									 PinitMag*sintheta*sin(phi),
									 PinitMag*costheta);

		

		// select ROI
		if(! (Ex>Ex_min && Ex<Ex_max) ) continue;

		// DOIT 
		//deex->DoDeex(Z,N,Ex); // this can also work w/ zero momentum
		deex->DoDeex(Z,N,Ex,Pinit);

		// scoling
		vector<Particle> particle = deex->GetParticleVector();
		eventID = deex->GetEventID();
		PinitX   = Pinit.X();
		PinitY   = Pinit.Y();
		PinitZ   = Pinit.Z();
		
		size=particle.size();
		for(int i=0;i<size;i++){
			Particle p = particle.at(i);
			PDG[i]=p._PDG;
			mass[i]=p._mass;
			totalE[i]=p.totalE();
			kE[i]=p.kE();
			PMag[i]=p._momentum.Mag();
			PX[i]=p._momentum.X();
			PY[i]=p._momentum.Y();
			PZ[i]=p._momentum.Z();
		}
		tree->Fill();
	}

	gStyle->SetTextSize(0.08);
	gStyle->SetTitleSize(0.045);
	gStyle->SetTitleXSize(0.045);
	gStyle->SetTitleYSize(0.045);
	gStyle->SetTitleYOffset(0.95);

	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	c->Divide(2,2);
	c->cd(1);
	h_sf_Ex_random->GetXaxis()->SetRangeUser(8,47);
	h_sf_Ex_random->SetStats(0);
	h_sf_Ex_random->SetMinimum(0);
	h_sf_Ex_random->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_sf_Ex_random->Draw("HIST");
	TLine* l_Ex_min = new TLine(Ex_min,0,Ex_min,h_sf_Ex_random->GetMaximum()*1.05);
	l_Ex_min->SetLineStyle(2);
	l_Ex_min->SetLineColor(kRed);
	l_Ex_min->Draw("same");
	TLine* l_Ex_max = new TLine(Ex_max,0,Ex_max,h_sf_Ex_random->GetMaximum()*1.05);
	l_Ex_max->SetLineStyle(2);
	l_Ex_max->SetLineColor(kRed);
	l_Ex_max->Draw("same");
	//
	c->cd(2);
	gPad->SetLogz();
	h_sf_random->SetStats(0);
	h_sf_random->SetMinimum(1);
	h_sf_random->Draw("colz");
	h_sf_random->SetStats(0);
	c->cd(3);
	h_sf_E_random->SetStats(0);
	h_sf_E_random->SetMinimum(0);
	h_sf_E_random->GetXaxis()->SetRangeUser(0,100);
	h_sf_E_random->Draw("HIST");
	c->cd(4);
	h_sf_p_random->SetStats(0);
	h_sf_p_random->SetMinimum(0);
	h_sf_p_random->Draw("HIST");
	c->Print("fig_sim/fig_sf.pdf");


	outf->cd();
	h_sf_Ex_random->Write();
	h_sf_E_random->Write();
	h_sf_p_random->Write();
	h_sf_random->Write();
	tree->Write();
	outf->Close();
	delete outf;



	return 0;
}
