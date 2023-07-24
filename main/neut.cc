#include <iostream>
#include <string>

#include "neutvtx.h"
#include "neutvect.h"
#include "neutpart.h"

#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "THStack.h"

#include "NucleusTable.hh"
#include "Deexcitation.hh"

using namespace std;
int main(int argc, char* argv[]){
	//------------------//
	string prefix = "neut_1GeV_numu_CCQE_C_MDLQE422";
	if(argc==2) prefix = argv[1];
	int seed=1;
	//-----------------.//

	int Zt, Nt;
	double S;
	ostringstream os;

	// prepare deex tools
	Deexcitation* deex = new Deexcitation(2, 1);
	deex->SetSeed(seed);
	deex->SetVerbose(1);
  NucleusTable* nucleus_table = deex->GetNucleusTablePtr();
	if(prefix.find("_C_")){
		Zt=6;
		Nt=6;
		os << getenv("TALYS_WORK_TABLES") << "/sf/pke12_tot.root";
		S = nucleus_table->GetNucleusPtr("12C")->S[2];
	}else if(prefix.find("_O_")){
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

	//--- input neut root
	os.str("");
	os << "output_neut/" << prefix.c_str() << ".root";
  TFile* rootf = new TFile(os.str().c_str(),"READ");
  
  TTree  *tree = (TTree *)(rootf->Get("neuttree"));
  NeutVtx *nvtx = new NeutVtx();
  NeutVect *nvect = new NeutVect();
  tree->SetBranchAddress("vertexbranch",&nvtx);
  tree->SetBranchAddress("vectorbranch",&nvect);
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


  int nevents = tree->GetEntries()/50;

  for ( int j = 0 ; j < nevents ; j++ ){
		tree->GetEntry(j);
		cout << "---------------------------------------------" << "\n";
		cout << "Event #        :" << nvect->EventNo << "\n";
		cout << "Target A       :" << nvect->TargetA << "\n";
		cout << "Target Z       :" << nvect->TargetA << "\n";
		cout << "VNuclIni       :" << nvect->VNuclIni << "\n";
		cout << "VNuclFin       :" << nvect->VNuclFin << "\n";
		cout << "PF Surface     :" << nvect->PFSurf   << "\n";
		cout << "PF Maximum     :" << nvect->PFMax    << "\n";
		cout << "Flux ID        :" << nvect->FluxID   << "\n";

		cout << "Intr. mode     :" << nvect->Mode   << "\n";

		double Ev=-1, Evv=-1, Enucc=-1;
		double Massnuc=-1;
		int NuPID=0;
		int nmulti=0;
		int Z=Zt,N=Nt;
		for ( int i = 0 ; i < nvect->Npart() ; i++ ){
			NeutPart* part = nvect->PartInfo(i);
			int fPID     = part->fPID;
			double fMass = part->fMass;
			int fStatus  = part->fStatus;
			bool fIsAlive = part->fIsAlive;

			double kE = part->fP.E()-fMass;

			// Init particles
			if(fStatus==-1 && fIsAlive==0){
				if(abs(fPID)==14){ //init nu
					Ev=kE;
					NuPID=fPID;
				}else{// target nuc
					Massnuc=fMass; 
				}
			}else if(fIsAlive==0){ // intermediate
				Enucc=part->fP.E();
			}else if(fIsAlive==1){ // postFSI
				if(abs(fPID)==13) Evv=part->fP.E(); // charged lepton
				else if(Enucc<0 && ( (NuPID>0 && fPID==2212) || (NuPID<0 && fPID==2112) )) Enucc = part->fP.E();// fsi target nuc
				if(fPID==2212) Z--;
				else if(fPID==2112){
					N--;
					nmulti++;
				}
			}
			cout << "i=" << i << "\n";
			cout << "Vertex         =" << nvect->VertexID(i) << "\n";
			cout << "Parent Index   =" << nvect->ParentIdx(i) << "\n";

			cout << "Particle Code  = " << (nvect->PartInfo(i))->fPID   << "\n";
			cout << "Particle Mass  = " << (nvect->PartInfo(i))->fMass   << "\n";
			cout << "Particle Mom.  =(" << (nvect->PartInfo(i))->fP.Px() << "," 
					 << (nvect->PartInfo(i))->fP.Py() << "," 
				 << (nvect->PartInfo(i))->fP.Pz() << "," 
				 << (nvect->PartInfo(i))->fP.E()  << ")" 
				 << "\n";
			cout << "Particle Flag  = " << (nvect->PartInfo(i))->fIsAlive << "\n";
			cout << "Particle Stat. = " << (nvect->PartInfo(i))->fStatus  << "\n";
			cout << "Particle Pos(1)=(" << (nvect->PartInfo(i))->fPosIni.X() << "," 
					 << (nvect->PartInfo(i))->fPosIni.Y() << "," 
				 << (nvect->PartInfo(i))->fPosIni.Z() << "," 
				 << (nvect->PartInfo(i))->fPosIni.T()  << ")" 
				 << "\n";
			cout << "Particle Pos(2)=(" << (nvect->PartInfo(i))->fPosFin.Px() << "," 
					 << (nvect->PartInfo(i))->fPosFin.Y() << "," 
				 << (nvect->PartInfo(i))->fPosFin.Z() << "," 
				 << (nvect->PartInfo(i))->fPosFin.T()  << ")" 
				 << "\n";
		}//end of part loop

		h_nmulti_postFSI->Fill(nmulti);
		double MissE = Ev-Evv-Enucc+Massnuc;
		h_MissE->Fill(MissE);

		double Ex = MissE-S;
		cout << MissE << "   "  << S <<  "  " << Ex << endl;
		h_Ex[0]->Fill(Ex);
		if(nvect->Npart()==4){
			deex->DoDeex(Zt,Nt,Z,N,0,Ex,TVector3(0,0,0));
			int shell = deex->GetShell();
			h_Ex[shell]->Fill(Ex);

			// add these elements to the neut output
			vector<Particle>* particle = deex->GetParticleVector();
			int size=particle->size();
			for(int i=0;i<size;i++){
				Particle p = particle->at(i);
				if(p._PDG==2112) nmulti++;
			}
			h_nmulti_postdeex->Fill(nmulti);
		}
	}





	TCanvas* c_nmulti_postFSI = new TCanvas("c_nmulti_postFSI","c_nmulti_postFSI",0,0,800,600);
	h_nmulti_postFSI->GetXaxis()->SetTitle("Neutron multiplicity");
	h_nmulti_postFSI->GetYaxis()->SetTitle("Events/bin");
	h_nmulti_postFSI->GetYaxis()->SetMaxDigits(2);
	h_nmulti_postFSI->Draw("HIST");
	h_nmulti_postdeex->SetLineColor(kRed);
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
	os << "fig_neut/fig_nmulti_postFSI_" << prefix.c_str() << ".pdf";
	c_nmulti_postFSI->Print(os.str().c_str());


	TCanvas* c_MissE = new TCanvas("c_MissE","c_MissE",0,0,800,600);
	h_MissE->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_MissE->GetYaxis()->SetTitle("Events/bin");
	h_MissE->GetYaxis()->SetMaxDigits(2);
	h_MissE->Draw("HIST");
	os.str("");
	os << "fig_neut/fig_MissE_" << prefix.c_str() << ".pdf";
	c_MissE->Print(os.str().c_str());

	TCanvas* c_Ex = new TCanvas("c_Ex","c_Ex",0,0,800,600);
	h_Ex[0]->GetXaxis()->SetTitle("Missing energy (MeV)");
	h_Ex[0]->GetYaxis()->SetTitle("Events/bin");
	h_Ex[0]->GetYaxis()->SetMaxDigits(2);
	h_Ex[0]->Draw("HIST");
	THStack* h_s_Ex = new THStack("h_s_Ex","");
	const int color[4]={1,600-7,632-7,920};
	for(int i=1;i<4;i++){
		h_Ex[i]->SetFillColor(color[i]);
		h_s_Ex->Add(h_Ex[i]);
	}
	h_s_Ex->Draw("same");
	gPad->RedrawAxis();
	os.str("");
	os << "fig_neut/fig_Ex_" << prefix.c_str() << ".pdf";
	c_Ex->Print(os.str().c_str());


	os.str("");
	os << "output_neut/histogram_deex_" << prefix.c_str() << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	h_nmulti_postFSI->Write();
	h_MissE->Write();
	outf->Close();
	delete outf;

	rootf->Close();
	delete rootf;

	return 0;
}
