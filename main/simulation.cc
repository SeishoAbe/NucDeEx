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
#include <TArrow.h>
#include <TText.h>
#include <TPaveText.h>


using namespace std;

int main(int argc, char* argv[]){
	if(argc<=1){
		cerr << "Input: " << argv[0] << " [Target nucleus] [Random Seed (optional)]" << endl;
		return 0;
	}
	int seed=1; // default: 1
	if(argc==3) seed = atoi(argv[2]);

	// ---- FIXME --- // 
	const int numofevent=1e5; // to be generated
	const double Ex_min =-1;
	const double Ex_max =-1;
		// negative -> not applied
	const bool flag_fig=0;
	// -------------- //

	ostringstream os,os_remove_g;

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
	os << "sim_out/" << argv[1] << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	TTree* tree = new TTree("tree",""); // LAB Freme, MeV
	int eventID=0, size;
	double MissE, Ex, S;
	double PinitMag, PinitX,PinitY,PinitZ;
	//
	int PDG[bins];
	double mass[bins];
	double totalE[bins],kE[bins];
	double PMag[bins], PX[bins],PY[bins],PZ[bins];
	string decay, decay_remove_g;
	tree->Branch("eventID",&eventID,"eventID/I");
	tree->Branch("decay",&decay);
	tree->Branch("decay_remove_g",&decay_remove_g);
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
	os.str("");	
	if(Z+N==11){ // 12C
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke12_tot.root";
		S = nucleus_table->GetNucleusPtr("12C")->S[2];
	}else if(Z+N==15){
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke16.root";
		S = nucleus_table->GetNucleusPtr("16O")->S[2];
	}
	cout << "S = " << S << endl;
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	TH2D* h_sf_int = (TH2D*) rootf->Get("h_sf_int");
	cout << gRandom->GetSeed() << endl;
	gRandom->SetSeed(seed); // for TH2 GetRandom2



	TCanvas* c_detail = new TCanvas("c_detail","",0,0,800,800);
	string pdfname = "fig_sim/fig_detail.pdf";
	if(flag_fig){
		c_detail->Print( (pdfname + (string)"[").c_str() );
		c_detail->Update();
		c_detail->Clear();
	}
	
	while(eventID<numofevent){
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
		if(Ex<0) continue;
		if(Ex_min>0 && Ex<Ex_min) continue;
		if(Ex_max>0 && Ex>Ex_max) continue;

		// DOIT 
		//deex->DoDeex(Z,N,Ex); // this can also work w/ zero momentum
		deex->DoDeex(Z,N,Ex,Pinit);

		// scoling
		vector<Particle> particle = deex->GetParticleVector();
		PinitX   = Pinit.X();
		PinitY   = Pinit.Y();
		PinitZ   = Pinit.Z();
		
		size=particle.size();
		os.str("");
		os_remove_g.str("");
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
			string pname;
			if(p._name.length()>4) pname = p._name.substr(0,1);
			else pname = p._name;
			os << pname.c_str();
			if(i!=0 && pname == "g") continue; // remove g if it is second
			os_remove_g << pname.c_str();
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
				for(int p=0;p<num_particle;p++){
					if(PDG[i] == PDG_particle[p]) color[i]=color_root[p];
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
				lXY[i]->Draw();
			}
			//
			c_detail->cd(2);
			TH1F* waku2= gPad->DrawFrame(-200,-200,200,200);
			waku2->SetTitle("YZ plane");
			waku2->GetXaxis()->SetTitle("Py (MeV)");
			waku2->GetYaxis()->SetTitle("Pz (MeV)");
			linitYZ->Draw();
			for(int i=0;i<size;i++){
				lYZ[i]->Draw();
			}
			//
			c_detail->cd(3);
			TH1F* waku3= gPad->DrawFrame(-200,-200,200,200);
			waku3->SetTitle("XZ plane");
			waku3->GetXaxis()->SetTitle("Px (MeV)");
			waku3->GetYaxis()->SetTitle("Pz (MeV)");
			linitXZ->Draw();
			for(int i=0;i<size;i++){
				lXZ[i]->Draw();
			}
			//
			c_detail->cd(4);
			TPaveText* t = new TPaveText(0.1,0.1,0.9,0.9);
			os.str("");
			os << "eventID = " << eventID;
			t->AddText(os.str().c_str());
			os.str("");
			os << argv[1] << ": Ex = " << fixed << setprecision(2) << Ex;
			t->AddText(os.str().c_str());
			os.str("");
			os << "P=(" << setprecision(1) << PinitX
				 << ", " << setprecision(1) << PinitY
				 << ", " << setprecision(1) << PinitZ << ")";
			t->AddText(os.str().c_str());
			//
			t->AddText("--- Decay products ---");
			for(int i=0;i<size;i++){
				os.str("");
				os << decay.c_str();
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
	if(flag_fig)c_detail->Print( (pdfname + (string)"]").c_str() );
	delete c_detail;


	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	c->Divide(2,2);
	c->cd(1);
	h_sf_Ex_random->GetXaxis()->SetRangeUser(-10,150);
	h_sf_Ex_random->SetStats(0);
	h_sf_Ex_random->SetMinimum(0);
	h_sf_Ex_random->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_sf_Ex_random->GetYaxis()->SetTitle("Events/bin");
	h_sf_Ex_random->Draw("HIST");
	TLine* l_Ex_min = new TLine(Ex_min,0,Ex_min,h_sf_Ex_random->GetMaximum()*1.05);
	l_Ex_min->SetLineStyle(2);
	l_Ex_min->SetLineColor(kRed);
	if(Ex_min>0) l_Ex_min->Draw("same");
	TLine* l_Ex_max = new TLine(Ex_max,0,Ex_max,h_sf_Ex_random->GetMaximum()*1.05);
	l_Ex_max->SetLineStyle(2);
	l_Ex_max->SetLineColor(kRed);
	if(Ex_max>0) l_Ex_max->Draw("same");
	//
	c->cd(2);
	gPad->SetLogz();
	h_sf_random->SetStats(0);
	h_sf_random->SetMinimum(1);
	h_sf_random->GetXaxis()->SetTitle("Momentum (MeV)");
	h_sf_random->GetYaxis()->SetTitle("Excitation energy (MeV)");
	h_sf_random->Draw("colz");
	h_sf_random->SetStats(0);
	c->cd(3);
	h_sf_E_random->SetStats(0);
	h_sf_E_random->SetMinimum(0);
	h_sf_E_random->GetXaxis()->SetRangeUser(0,150);
	h_sf_E_random->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_sf_E_random->GetYaxis()->SetTitle("Events/bin");
	h_sf_E_random->Draw("HIST");
	TLine* l_S = new TLine(S,0,S,h_sf_E_random->GetMaximum()*1.05);
	l_S->SetLineStyle(2);
	l_S->SetLineColor(kRed);
	l_S->Draw("same");
	c->cd(4);
	h_sf_p_random->SetStats(0);
	h_sf_p_random->SetMinimum(0);
	h_sf_p_random->GetXaxis()->SetTitle("Momentum (MeV)");
	h_sf_p_random->GetYaxis()->SetTitle("Events/bin");
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
