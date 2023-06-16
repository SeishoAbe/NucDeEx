#include "include/consts.hh"
#include <TStyle.h>
#include <map>

using namespace std;
int plot_simulation(){
	ostringstream os;
	TFile* rootf = new TFile("sim_out/11B.root","READ");
	TTree* tree = (TTree*) rootf->Get("tree");
	int eventID, size;
	double MissE, Ex, S;
	double PinitMag, PinitX,PinitY,PinitZ;
	//
	int PDG[bins];
	double mass[bins];
	double totalE[bins],kE[bins];
	double PMag[bins], PX[bins],PY[bins],PZ[bins];
	string* decay=0;
	tree->SetBranchAddress("eventID",&eventID);
	//tree->SetBranchAddress("decay",&decay);
	tree->SetBranchAddress("decay_remove_g",&decay);
	tree->SetBranchAddress("MissE",&MissE);
	tree->SetBranchAddress("S",&S);
	tree->SetBranchAddress("Ex",&Ex);
	tree->SetBranchAddress("PinitMag",&PinitMag);
	tree->SetBranchAddress("PinitX",&PinitX);
	tree->SetBranchAddress("PinitY",&PinitY);
	tree->SetBranchAddress("PinitZ",&PinitZ);
	tree->SetBranchAddress("size",&size);
	tree->SetBranchAddress("PDG",&PDG);
	tree->SetBranchAddress("mass",&mass);
	tree->SetBranchAddress("totalE",&totalE);
	tree->SetBranchAddress("kE",&kE);
	tree->SetBranchAddress("PMag",&PMag);
	tree->SetBranchAddress("PX",&PX);
	tree->SetBranchAddress("PY",&PY);
	tree->SetBranchAddress("PZ",&PZ);

	// --- Draw --- // 
	gStyle->SetTextSize(0.08);
	gStyle->SetTitleSize(0.045);
	gStyle->SetTitleXSize(0.045);
	gStyle->SetTitleYSize(0.045);
	gStyle->SetTitleYOffset(0.95);
	
	TH1D* h_nmulti = new TH1D("h_nmulti","",10,-0.5,9.5);
	TH1D* h_kE[num_particle];
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "h_kE_" << p;
		h_kE[p]= new TH1D(os.str().c_str(),"",40,0,20);
	}

	int max_size=0;
	int numofevent=tree->GetEntries();
	cout << "numofevent=" << numofevent << endl;
	map<string, double> br;
	map<string, double> :: iterator itr;
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
		if(max_size<size) max_size=size;
		int nmulti=0;
		for(int b=0;b<size;b++){
			if(PDG[b]==2112) nmulti++;
			int p=-1;
			for(int par=0;par<num_particle;par++){
				if(PDG_particle[par]==PDG[b]){
					p=par;
					break;
				}
			}
			if(p>=0) h_kE[p]->Fill(kE[b]);
		}
		itr = br.find(decay->c_str());
		if(itr != end(br) ) {
			itr->second += 1;
		}else{
			cout << decay->c_str() << endl;
			br.insert(make_pair(decay->c_str(),1));
		}
		h_nmulti->Fill(nmulti);
	}
	cout << "max_size=" << max_size << endl;
	for(itr = br.begin();itr!=br.end();itr++){
		itr->second /= tree->GetEntries()/100; // %
		cout << setw(10) << itr->first << "  " << itr->second << endl;
	}

	// Scale
	h_nmulti->Scale(1./h_nmulti->GetEntries());
	for(int p=0;p<num_particle;p++){
		h_kE[p]->Scale(1./h_kE[p]->GetEntries());
	}


	TCanvas* c_nmulti = new TCanvas("c_nmulti","",0,0,800,600);
	h_nmulti->GetXaxis()->SetNdivisions(10);
	h_nmulti->GetXaxis()->SetTitle("Neutron multiplicity");
	h_nmulti->GetYaxis()->SetTitle("Branching ratio");
	h_nmulti->SetStats(0);
	h_nmulti->Draw("HIST");
	os.str("");
	os << "Mean = " << setprecision(3) << h_nmulti->GetMean();
	TText* t_nmulti = new TText(5,h_nmulti->GetMaximum()*0.5,os.str().c_str());
	t_nmulti->Draw("same");
	c_nmulti->Print("fig_sim/fig_nmulti.pdf");
	
	TCanvas* c_kE = new TCanvas("c_kE","",0,0,1200,600);
	c_kE->Divide(4,2);
	for(int p=0;p<num_particle;p++){
		c_kE->cd(p+1);
		h_kE[p]->SetTitle(particle_name[p].c_str());
		h_kE[p]->SetStats(0);
		h_kE[p]->GetXaxis()->SetTitle("Kinetic energy (MeV)");
		h_kE[p]->GetYaxis()->SetTitle("A.U.");
		h_kE[p]->Draw("HIST");
		os.str("");
		os << "Prob = " << fixed << (double)h_kE[p]->GetEntries()/tree->GetEntries();
		TText* t_kE = new TText(7,h_kE[p]->GetMaximum()*0.7,os.str().c_str());
		t_kE->Draw("same");
	}
	c_kE->cd(8);
	TPaveText* t = new TPaveText(0.1,0.1,0.9,0.9);
	os.str("");
	os << "# generated events = " << numofevent;
	t->AddText(os.str().c_str());
	t->Draw("same");
	c_kE->Print("fig_sim/fig_kE.pdf");


	return 0;
}
