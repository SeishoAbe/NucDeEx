#include <TStyle.h>
#include <algorithm>
#include <map>

#include "../include/consts.hh"
#include "../include/NucleusTable.hh"
#include "../include/Nucleus.hh"

R__LOAD_LIBRARY(lib/libTALYStool);

using namespace std;

int plot_simulation(){
	// ---- FIXME ---- //
	string target = "11B";
	const double Ex_min =16;
	const double Ex_max =35;
		// negative -> not applied
	const int ldmodel=1;
	const bool parity_optmodall=1;
	const bool flag_decay=1;
		// 0 -> use string w/ "g" (gamma)
		// 1- > use string w/o "g" (gamma) <- use this 
	// --------------- //

	ostringstream os;
	os.str("");
	os << "sim_out/" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << ".root";
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	cout << os.str().c_str() << endl;
	TTree* tree = (TTree*) rootf->Get("tree");
	int eventID, size;
	double MissE, Ex, S;
	double PinitMag, PinitX,PinitY,PinitZ;
	//
	int PDG[bins];
	double mass[bins];
	double totalE[bins],kE[bins];
	double PMag[bins], PX[bins],PY[bins],PZ[bins];
	bool flag[bins];
	double Ex_daughter[bins];
	string* decay=0;
	tree->SetBranchAddress("eventID",&eventID);
	if(flag_decay) tree->SetBranchAddress("decay_remove_g",&decay);
	else           tree->SetBranchAddress("decay",&decay);
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
	tree->SetBranchAddress("flag",&flag);
	tree->SetBranchAddress("Ex_daughter",&Ex_daughter);

	// --- Draw --- // 
	gStyle->SetTextFont(132);
	gStyle->SetTextSize(0.08);
	gStyle->SetTitleSize(0.05,"XYZ");
	gStyle->SetTitleFont(132,"XYZ");
	gStyle->SetLabelSize(0.05,"XYZ");
	gStyle->SetLabelFont(132,"XYZ");
	gStyle->SetLegendFont(132);
	gStyle->SetTitleYOffset(0.95);
	
	TH1D* h_Ex = new TH1D("h_Ex","",500,-100,400);
	TH1D* h_nmulti = new TH1D("h_nmulti","",10,-0.5,9.5);
	TH1D* h_kE[num_particle];
	TH1D* h_Ex_particle[num_particle]; // w/o experimental energy th
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "h_kE_" << p;
		if(p==0) h_kE[p]= new TH1D(os.str().c_str(),"",200,0,20);
		else     h_kE[p]= new TH1D(os.str().c_str(),"",40,0,20);
	}

	int max_size=0;
	int numofevent=0;
	map<string, double> br;
	map<string, double> :: iterator itr;
	// for comparison w/ Panin
	double rbr[num_particle]={0}; // two-body
	double rbr_2b[num_particle]={0}; // two-body
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);

		// --- Ex selection --- //
		h_Ex->Fill(Ex);
		if(Ex_min>0 && Ex<Ex_min) continue;
		if(Ex_max>0 && Ex>Ex_max) continue;
		numofevent++;

		if(max_size<size) max_size=size;

		int multi[num_particle]={0};
		int particle_counter=0;
		int first_p=-1;
		for(int b=0;b<size;b++){
			int p=-1;
			for(int par=0;par<num_particle;par++){
				if(PDG_particle[par]==PDG[b]){
					p=par;
					break;
				}
			}
			if(p>=0){ // this is particle info (not daughter nuclei)
				multi[p]++;
				h_kE[p]->Fill(kE[b]);
				h_Ex_particle[p]->Fill(Ex);
				if(first_p<0 && p>=1) first_p=p; // except gamma
				if(p>=1) particle_counter++;
			}
		}

		// --- Fill multiplicity histogram
		h_nmulti->Fill(multi[1]); // w/o energy th
		
		// --- Nominal rbr 
		if(first_p>=0){
			rbr[first_p]++;
			if(particle_counter==1) rbr_2b[first_p]++;
		}

		//--- save map--//
		itr = br.find(decay->c_str());
		if(itr != end(br) ) {
			itr->second += 1;
		}else{
			//cout << decay->c_str() << endl;
			br.insert(make_pair(decay->c_str(),1));
		}
	}
	cout << "max_size=" << max_size << endl;
	cout << "numofevent=" << numofevent << endl;

	cout << "### without threshold " << endl;
	double rbr_sum=0, rbr_2b_sum=0;
	for(int par=0;par<num_particle;par++){
		rbr[par] = rbr[par]/numofevent*100; // %
		rbr_2b[par] = rbr_2b[par]/numofevent*100; // %
		cout << setw(10) << particle_name[par].c_str() << "  " 
				 << setw(10) << fixed << setprecision(1) << rbr[par] << "  " << rbr_2b[par] << endl;
		rbr_sum+=rbr[par];
		rbr_2b_sum+=rbr_2b[par];
	}
	cout << setw(10) << "sum" << "  " << setw(10) << rbr_sum << "  " << rbr_2b_sum << endl;

	
	// output map as .txt
	os.str("");
	if(Ex_min>0) os << "_Exmin" << Ex_min;
	if(Ex_max>0) os	<< "_Exmax" << Ex_max;
	string suffix = os.str();
	os.str("");
	os << "sim_out/" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << suffix.c_str() << "_raw.txt";
	ofstream ofs_raw (os.str().c_str());
	os.str("");
	os << "sim_out/" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << suffix.c_str() << ".txt";
	ofstream ofs (os.str().c_str());
	const double br_th = 0.5;
	ofs << "# br threshold written in this table: " << br_th << endl;
	double br_sum=0; // for check
	double br_others=0;
	for(itr = br.begin();itr!=br.end();itr++){
		double br_s = itr->second/numofevent*100;
		br_sum += br_s;
		ofs_raw  << setw(10) << itr->first << "   " 
				 << setw(10) << fixed << setprecision(2) << br_s << " (%)" << endl;
		if(br_s>br_th){
			ofs  << setw(10) << itr->first << "   " 
					 << setw(10) << fixed << setprecision(2) << br_s << " (%)" << endl;
		}else br_others+=br_s;
	}
	cout << "(check) br_sum = " << br_sum << " (%)" << endl;
	ofs_raw.close();
	ofs  << setw(10) << "others" << "   " 
			 << setw(10) << fixed << setprecision(2) << br_others << " (%)" << endl;
	ofs.close();

	// --- Scale histograms 
	h_Ex->Scale(1./h_Ex->GetEntries());
	h_nmulti->Scale(1./h_nmulti->GetEntries());
	for(int p=0;p<num_particle;p++){
		h_kE[p]->Scale(1./h_kE[p]->GetEntries());
		h_Ex_particle[p]->Scale(2.0/h_Ex->GetEntries());
	}

	TCanvas* c_Ex = new TCanvas("c_Ex","",0,0,800,600);
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_Ex" << suffix.c_str() << ".pdf";
	string pdfname=os.str();
	c_Ex->Print( (pdfname+"[").c_str());
	c_Ex->Update();
	c_Ex->Clear();
	//
	h_Ex->GetXaxis()->SetRangeUser(0,100);
	h_Ex->SetStats(0);
	h_Ex->SetMinimum(0);
	h_Ex->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_Ex->GetYaxis()->SetTitle("A.U.");
	h_Ex->Draw("HIST");
	TLine* l_1_Ex_min = new TLine(Ex_min,0,Ex_min,h_Ex->GetMaximum()*1.05);
	l_1_Ex_min->SetLineWidth(2);
	l_1_Ex_min->SetLineStyle(2);
	l_1_Ex_min->SetLineColor(kRed);
	if(Ex_min>8) l_1_Ex_min->Draw("same");
	TLine* l_1_Ex_max = new TLine(Ex_max,0,Ex_max,h_Ex->GetMaximum()*1.05);
	l_1_Ex_max->SetLineWidth(2);
	l_1_Ex_max->SetLineStyle(2);
	l_1_Ex_max->SetLineColor(kRed);
	if(Ex_max>8) l_1_Ex_max->Draw("same");
	c_Ex->Print(pdfname.c_str());
	c_Ex->Update();
	c_Ex->Clear();
	//
	h_Ex->SetTitle("w/o energy th. (double count is included)");
	h_Ex->GetXaxis()->SetRangeUser(8,47);
	h_Ex->Draw("HIST");
	for(int p=1;p<num_particle;p++){
		h_Ex_particle[p]->Draw("HISTsame");
	}
	TLine* l_Ex_min = new TLine(Ex_min,0,Ex_min,h_Ex->GetMaximum()*1.05);
	l_Ex_min->SetLineWidth(2);
	l_Ex_min->SetLineStyle(2);
	l_Ex_min->SetLineColor(kRed);
	if(Ex_min>8 && Ex_min<47) l_Ex_min->Draw("same");
	TLine* l_Ex_max = new TLine(Ex_max,0,Ex_max,h_Ex->GetMaximum()*1.05);
	l_Ex_max->SetLineWidth(2);
	l_Ex_max->SetLineStyle(2);
	l_Ex_max->SetLineColor(kRed);
	if(Ex_max>8 && Ex_max<47) l_Ex_max->Draw("same");
	c_Ex->Print(pdfname.c_str());
	c_Ex->Update();
	c_Ex->Clear();
	//
	THStack* h_s_Ex = new THStack("h_s_Ex","");
	h_Ex->Draw("HIST");
	for(int p=2;p<num_particle;p++){
		h_Ex_particle[p]->SetFillColor(color_root[p]);
		h_s_Ex->Add(h_Ex_particle[p]);
	}
	h_s_Ex->Draw("HISTsame");
	c_Ex->Print((pdfname+"]").c_str());
	c_Ex->Update();
	c_Ex->Clear();


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
	//
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_nmulti" << suffix.c_str() << ".pdf";
	c_nmulti->Print(os.str().c_str());

	
	TCanvas* c_kE = new TCanvas("c_kE","",0,0,1200,600);
	c_kE->Divide(4,2);
	for(int p=0;p<num_particle;p++){
		c_kE->cd(p+1);
		h_kE[p]->SetTitle(particle_name[p].c_str());
		h_kE[p]->SetStats(0);
		h_kE[p]->GetXaxis()->SetTitle("Kinetic energy (MeV)");
		h_kE[p]->GetYaxis()->SetTitle("A.U.");
		h_kE[p]->SetLineColor(color_root[p]);
		if(p!=0) h_kE[p]->SetLineWidth(2);
		h_kE[p]->Draw("HIST");
		os.str("");
		os << "Prob = " << fixed << (double)h_kE[p]->GetEntries()/numofevent;
		TText* t_kE = new TText(7,h_kE[p]->GetMaximum()*0.7,os.str().c_str());
		t_kE->Draw("same");
	}
	c_kE->cd(8);
	TPaveText* t = new TPaveText(0.1,0.1,0.9,0.9);
	os.str("");
	os << "# generated events = " << numofevent;
	t->AddText(os.str().c_str());
	t->Draw("same");
	//
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_kE" << suffix.c_str() << ".pdf";
	c_kE->Print(os.str().c_str());


	return 0;
}
