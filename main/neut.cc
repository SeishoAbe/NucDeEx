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

#include "neutvtx.h"
#include "neutvect.h"
#include "neutpart.h"

#include "NucleusTable.hh"
#include "Deexcitation.hh"

using namespace std;
int main(int argc, char* argv[]){
	//------------------//
	//string prefix = "neut_1GeV_numu_CCQE_C_MDLQE422";
	//string prefix = "neut_1GeV_numub_CCQE_C_MDLQE422";
	//string prefix = "neut_1GeV_numu_CCQE_O_MDLQE422_NUCDEXITE0";
	//string prefix = "neut_1GeV_numub_CCQE_O_MDLQE422_NUCDEXITE0";
	string prefix = "neut_1GeV_numu_NCQE_O_MDLQE422_NUCDEXITE0";
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
	bool flag_O=0;
	if(prefix.find("_C_")!=string::npos){
		Zt=6;
		Nt=6;
		os << getenv("NUCDEEX_TABLES") << "/sf/pke12_tot.root";
		S = nucleus_table->GetNucleusPtr("12C")->S[2];
	}else if(prefix.find("_O_")!=string::npos){
		Zt=8;
		Nt=8;
		os << getenv("NUCDEEX_TABLES") << "/sf/pke16.root";
		S = nucleus_table->GetNucleusPtr("16O")->S[2];
		flag_O=1;
	}
	cout << "Separation E = " << S << endl;
	TFile* root = new TFile(os.str().c_str(),"READ");
	TH2D* h_sf_int = (TH2D*) root->Get("h_sf_int");
	h_sf_int->SetDirectory(0);
	gRandom->SetSeed(seed); // for TH2 GetRandom2
	root->Close();
	delete root;

	bool flag_CC=1;//1: CC, 0: NC
	if(prefix.find("NCQE")!=string::npos) flag_CC=0;

	// --- input neut root
	os.str("");
	os << "output_neut/" << prefix.c_str() << ".root";
  TFile* rootf = new TFile(os.str().c_str(),"READ");
	cout << os.str().c_str() << endl;
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

	TH1D* h_Pinit = new TH1D("h_Pinit","",100,0,500);
	TH1D* h_nmulti_postFSI = new TH1D("h_nmulti_postFSI","",10,-0.5,9.5);
	TH1D* h_nmulti_postdeex = new TH1D("h_nmulti_postdeex","",10,-0.5,9.5);
	TH1D* h_MissE = new TH1D("h_MissE","",400,0,200);
	TH1D* h_Ex[4]; // single [0]: all
	for(int i=0;i<4;i++){
		os.str("");
		os << "h_Ex_" << i;
		h_Ex[i] = new TH1D(os.str().c_str(),"",500,-100,400);
	}
	TH1D* h_Ex_multi = new TH1D("h_Ex_multi","",500,-100,400); // multi nucleon hole
	TH2D* h_MissE_Pinit = new TH2D("h_MissE_Pinit","",100,0,500,400,0,200);

	// --- read root
  int nevents = tree->GetEntries();
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
		double Pinit=0;
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
					if(flag_CC){ // CC
						if(fPID==2112){// n->p
							Z++;
							N--;
						}else if(fPID==2212){//p->n
							Z--;
							N++;
						}
					}else{
						if(fPID==2112)N--;
						else if(fPID==2212) Z--;
					}
					Massnuc=fMass; 
					Pinit = sqrt( pow(nvect->PartInfo(i)->fP.Px(),2)
					             + pow(nvect->PartInfo(i)->fP.Py(),2)
					             + pow(nvect->PartInfo(i)->fP.Pz(),2) );
					h_Pinit->Fill(Pinit);
				}
			}else if(fIsAlive==0 && Enucc<0){ // intermediate
				Enucc=part->fP.E();
			}else if(fIsAlive==1){ // postFSI
				if(abs(fPID)==13 || abs(fPID)==14) Evv=part->fP.E(); // charged lepton
				else if(Enucc<0 && (fPID==2212 || fPID==2112))Enucc = part->fP.E();// fsi target nuc
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
		h_MissE_Pinit->Fill(Pinit,MissE);

		double Ex = MissE-S;
		//cout << MissE << "   "  << S <<  "  " << Ex << endl;
		h_Ex[0]->Fill(Ex);

		// --- DOIT --- //
		deex->DoDeex(Zt,Nt,Z,N,0,Ex,TVector3(0,0,0));

		// --- Scoring --- //
		int shell = deex->GetShell();
		h_Ex[shell]->Fill(Ex);
		vector<Particle>* particle = deex->GetParticleVector();
		int size=particle->size();
		for(int i=0;i<size;i++){
			Particle p = particle->at(i);
			if(p._PDG==2112) nmulti++;
		}
		h_nmulti_postdeex->Fill(nmulti);
  }


	// --- plot
	TCanvas* c_Pinit = new TCanvas("c_Pinit","c_Pinit",0,0,800,600);
	h_Pinit->GetXaxis()->SetTitle("Momentum of target nucleon (MeV)");
	h_Pinit->GetYaxis()->SetTitle("Events/bin");
	h_Pinit->GetYaxis()->SetMaxDigits(2);
	h_Pinit->Draw("HIST");
	os.str("");
	os << "fig_neut/fig_Pinit_" << prefix.c_str() << ".pdf";
	c_Pinit->Print(os.str().c_str());

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
	if(flag_O==1){
		os.str("");
		os << "Prob(p1/2)=" << fixed << setprecision(1) << (double)h_Ex[3]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
		TText* t_p12 = new TText(40,h_Ex[0]->GetMaximum()*0.7,os.str().c_str());
		t_p12->Draw("same");
	}
	//
	gPad->RedrawAxis();
	os.str("");
	os << "fig_neut/fig_Ex_" << prefix.c_str() << ".pdf";
	c_Ex->Print(os.str().c_str());

	TCanvas* c_MissE_Pinit = new TCanvas("c_MissE_Pinit","c_MissE_Pinit",0,0,800,600);
	gPad->SetLogz();
	h_MissE_Pinit->GetXaxis()->SetTitle("Momentum of target nucleon (MeV)");
	h_MissE_Pinit->GetYaxis()->SetTitle("Missing energy (MeV)");
	h_MissE_Pinit->GetYaxis()->SetRangeUser(0,100);
	h_MissE_Pinit->SetStats(0);
	h_MissE_Pinit->Draw("colz");
	os.str("");
	os << "fig_neut/fig_MissE_Pinit_" << prefix.c_str() << ".pdf";
	c_MissE_Pinit->Print(os.str().c_str());

	// --- save 
	os.str("");
	os << "output_neut/histogram_deex_" << prefix.c_str() << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	h_Pinit->Write();
	h_nmulti_postFSI->Write();
	h_nmulti_postdeex->Write();
	h_MissE->Write();
	for(int i=0;i<4;i++){
		h_Ex[i]->Write();
	}
	h_Ex_multi->Write();
	h_MissE_Pinit->Write();
	outf->Close();
	delete outf;

	rootf->Close();
	delete rootf;

	return 0;
}
