#include <string>
#include <iostream>

#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"

using namespace std;

const double mass_neutron=939.566;//MeV
const double mass_proton=938.272;//MeV

int main(int argc, char* argv[]){
	if(argc!=4){
		cerr << "Input: " << argv[0] << " [flavor] [target] [tune]" << endl;
		return 1;
	}

	// --- FIXME --- //
	const int flavor=atoi(argv[1]);
	string target = argv[2];
	string tune = argv[3];
	// ------------- //

	ostringstream os;
	os.str("");
	os << "CCQE." << flavor << "." << target.c_str() 
		 << "." << tune.c_str();
	string prefix = os.str();
	os.str("");
	os << "output/gntp." << prefix.c_str() << ".gst.root";
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	TTree* tree = (TTree*)rootf->Get("gst");
	const int kNPmax = 250;
	double Ev, pxv, pyv, pzv;
	int neu, tgt, hitnuc;
	double pxn, pyn, pzn;
	int ni;
	int pdgi[kNPmax];
	double Ei[kNPmax], pxi[kNPmax], pyi[kNPmax],pzi[kNPmax];
	int fspl;
	double El, pl, pxl, pyl, pzl;
	int nfn,nf;
	int pdgf[kNPmax];
	double Ef[kNPmax], pf[kNPmax], pxf[kNPmax], pyf[kNPmax], pzf[kNPmax];
	double sumKEf;
	// for nu
	tree->SetBranchAddress("Ev",&Ev);
	tree->SetBranchAddress("pxv",&pxv);
	tree->SetBranchAddress("pyv",&pyv);
	tree->SetBranchAddress("pzv",&pzv);
	tree->SetBranchAddress("neu",&neu);
	//fortarget&hitnuc
	tree->SetBranchAddress("tgt",&tgt);
	tree->SetBranchAddress("hitnuc",&hitnuc);
	tree->SetBranchAddress("pxn",&pxn);
	tree->SetBranchAddress("pyn",&pyn);
	tree->SetBranchAddress("pzn",&pzn);
	//forhadroninnucleus(='primary'hadronicsystem)
	tree->SetBranchAddress("ni",&ni);
	tree->SetBranchAddress("pdgi",pdgi);
	tree->SetBranchAddress("Ei",Ei);
	tree->SetBranchAddress("pxi",pxi);
	tree->SetBranchAddress("pyi",pyi);
	tree->SetBranchAddress("pzi",pzi);
	//forlepton
	tree->SetBranchAddress("fspl",&fspl);
	tree->SetBranchAddress("El",&El);
	tree->SetBranchAddress("pl",&pl);
	tree->SetBranchAddress("pxl",&pxl);
	tree->SetBranchAddress("pyl",&pyl);
	tree->SetBranchAddress("pzl",&pzl);
	//forfinalstate
	tree->SetBranchAddress("nfn",&nfn);
	tree->SetBranchAddress("nf",&nf);
	tree->SetBranchAddress("pdgf",pdgf);
	tree->SetBranchAddress("Ef",Ef);
	tree->SetBranchAddress("pf",pf);
	tree->SetBranchAddress("pxf",pxf);
	tree->SetBranchAddress("pyf",pyf);
	tree->SetBranchAddress("pzf",pzf);
	//forinteraction
	tree->SetBranchAddress("sumKEf",&sumKEf);

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
  TH1D* h_MissE = new TH1D("h_MissE","",400,0,200);


	//--- GeV2MeV --- //
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);

		h_Ev->Fill(Ev*1e3); 
		h_Pinit->Fill( sqrt(pxn*pxn+pyn*pyn+pzn*pzn)*1e3 );

		double Enucc=0;
		for(int k=0;k<ni;k++){
			if(pdgi[k]==2212 || pdgi[k]==2112) Enucc=Ei[k];
		}
		double massnuc=0;
		if(hitnuc==2112) massnuc = mass_neutron*1e-3; // MeV2GeV
		else if (hitnuc==2212) massnuc = mass_proton*1e-3; // MeV2GeV
		double MissE=(Ev-El-Enucc+massnuc)*1e3;
		h_MissE->Fill(MissE);
		//cout << MissE << "   " << Ev << "   " << El << "   " << Enucc
		//		 << "   " << massnuc << endl;


		int nmulti=0;
		for(int f=0;f<nf;f++){
			if(pdgf[f]==2112){//neutron
				h_kEn->Fill(Ef[f]*1e3-mass_neutron);
				h_lowkEn->Fill(Ef[f]*1e3-mass_neutron);
				nmulti++;
			}else if(pdgf[f]==2212){//proton
				h_kEp->Fill(Ef[f]*1e3-mass_proton);
				h_lowkEp->Fill(Ef[f]*1e3-mass_proton);
			}
			h_PID_postFSI->Fill(pdgf[f]);
		}
		h_nmulti_postFSI->Fill(nmulti);

		if(nmulti!=nfn){
			cerr << "something wrong happens.." << endl;
			abort();
		}
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
	h_MissE->GetXaxis()->SetRangeUser(0,100);
	h_MissE->Draw("HIST");
	os.str("");
	os << "fig/fig_MissE_" << prefix.c_str() << ".pdf";
	c_MissE->Print(os.str().c_str());



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
	outf->Close();





	return 0;
}
