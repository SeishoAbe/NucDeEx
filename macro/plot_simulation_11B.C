#include <TStyle.h>
#include <algorithm>
#include <map>

#include "../include/consts.hh"
#include "../include/NucleusTable.hh"
#include "../include/Nucleus.hh"

R__LOAD_LIBRARY(lib/libTALYStool);

using namespace std;

const double kE_th[num_particle]
	= {0,3.1,3.1,
		 4.0,4.6,0,4.5};

int plot_simulation_11B(){
	// ---- FIXME ---- //
	string target = "11B";
	const double Ex_min =16;
	const double Ex_max =35;
		// negative -> not applied
	const int ldmodel=2;
	const bool parity_optmodall=1;
	const bool flag_decay=1;
		// 0 -> use string w/ "g" (gamma)
		// 1- > use string w/o "g" (gamma) <- use this 
	// --------------- //

	NucleusTable* _nucleus_table = new NucleusTable();
	if(!_nucleus_table->ReadTables()){
		cerr << "Fatal Error" << endl;
		abort();
	}
	Nucleus* nuc_target = _nucleus_table->GetNucleusPtr(target.c_str());
	double SE[num_particle];
	string daughter[num_particle]={"11B","10B","10Be","9Be","8Be","8Li","7Li"};
	double criteria_3b[num_particle]={0};
	for(int p=0;p<num_particle;p++){
		SE[p] = nuc_target->S[p];
		Nucleus* nuc = _nucleus_table->GetNucleusPtr(daughter[p].c_str());
		criteria_3b[p] = min(5.0,(double)nuc->min_S());
		criteria_3b[p] += SE[p];
		cout << "criteria_3b = " << criteria_3b[p] << endl;
	}

	ostringstream os;
	os.str("");
	os << "sim_out/12C/" << target.c_str() << "_ldmodel" << ldmodel;
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
	gStyle->SetLegendTextSize(0.04);
	gStyle->SetTitleYOffset(0.95);
	
	TH1D* h_Ex = new TH1D("h_Ex","",500,-100,400);
	TH1D* h_nmulti = new TH1D("h_nmulti","",10,-0.5,9.5);
	TH1D* h_kE[num_particle];
	TH1D* h_Ex_particle[num_particle]; // w/o experimental energy th
	TH1D* h_Ex_particle_th[num_particle]; // w/ 
	TH2D* h_Ex_kE[num_particle]; 
	TH1D* h_SE[num_particle];
	for(int p=0;p<num_particle;p++){
		os.str("");
		os << "h_kE_" << p;
		if(p==0) h_kE[p]= new TH1D(os.str().c_str(),"",200,0,20);
		else     h_kE[p]= new TH1D(os.str().c_str(),"",40,0,20);
		//
		os.str("");
		os << "h_Ex_particle_" << p;
		h_Ex_particle[p] = new TH1D(os.str().c_str(),"",500,-100,400);
		h_Ex_particle[p]->SetLineColor(color_root[p]);
		//
		os.str("");
		os << "h_Ex_particle_th_" << p;
		h_Ex_particle_th[p] = new TH1D(os.str().c_str(),"",500,-100,400);
		h_Ex_particle_th[p]->SetLineColor(color_root[p]);
		//
		os.str("");
		os << "h_Ex_kE_" << p;
		h_Ex_kE[p] = new TH2D(os.str().c_str(),"",1000,-100,400,100,0,50);
		//
		os.str("");
		os << "h_SE_" << p;
		h_SE[p] = new TH1D(os.str().c_str(),"",50,0,50);
	}

	int max_size=0;
	int numofevent=0, numofdetected=0;
	map<string, double> br;
	map<string, double> :: iterator itr;
	// for comparison w/ Panin
	double rbr[num_particle]={0}; // two-body
	double rbr_2b[num_particle]={0}; // two-body
	// for comparison w/ Yosoi
	// double counter could be happen
	double rbr_th[num_particle]={0};  // w/ th
	double rbr_th_2b[num_particle]={0}; // w/ th two-body
	int count_particle_th[num_particle]={0};
	int count_particle_3b[num_particle]={0};
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);

		// --- Ex selection --- //
		h_Ex->Fill(Ex);
		if(Ex_min>0 && Ex<Ex_min) continue;
		if(Ex_max>0 && Ex>Ex_max) continue;
		numofevent++;

		if(max_size<size) max_size=size;

		int multi[num_particle]={0};
		bool flag_detected[num_particle]={0};
		int particle_counter=0,particle_counter_th=0; // w/o g
		int first_p=-1;
		double kE_r[num_particle]={0};
		double Ex_daughter_r=-1;
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
				//particle_counter++;
				// --- apply threshold 
				if(kE[b]<kE_th[p]) continue;
				// --- 
				kE_r[p]=kE[b];
				count_particle_th[p]++;
				flag_detected[p] = 1; // detected
				//if(p>=1) particle_counter_th++;
				particle_counter_th++;
			}else if(Ex_daughter_r<0 && first_p>0){ // daughter
				if(flag[b]==0) Ex_daughter_r=Ex_daughter[b];
				else Ex_daughter_r=0;
			}
		}

		// --- Fill multiplicity histogram
		h_nmulti->Fill(multi[1]); // w/o energy th
		
		// --- Nominal rbr 
		if(first_p>=0){
			rbr[first_p]++;
			if(particle_counter==1) rbr_2b[first_p]++;
		}
			
		// --- Detected rbr 
		bool detected=0;
		for(int par=1;par<num_particle;par++){
			if(!flag_detected[par]) continue;
			rbr_th[par]++;
			h_Ex_kE[par]->Fill(Ex,kE_r[par]);
			if(Ex_daughter_r>=0){
				h_SE[par]->Fill(Ex_daughter_r+SE[par]);
				if(Ex_daughter_r+SE[par]>criteria_3b[par]) count_particle_3b[par]++;
			}
			//if( (par!=1 &&particle_counter_th==1) || (par==1 && particle_counter==1) ){
			if(particle_counter_th==1){
				rbr_th_2b[par]++;
				h_Ex_particle_th[par]->Fill(Ex);
			}
			if(par>=2) detected=1; // neutrons did not affect detection
		}
		if(detected) numofdetected++;

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
	cout << "numofdetected=" << numofdetected << endl;

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

	cout << "### with threshold " << endl;
	double rbr_th_sum=0, rbr_th_2b_sum=0;
	for(int par=0;par<num_particle;par++){
		rbr_th[par] = rbr_th[par]/numofevent*100; // %
		rbr_th_2b[par] = rbr_th_2b[par]/numofevent*100; // %
		cout << setw(10) << particle_name[par].c_str() << "  " 
				 << setw(10) << rbr_th[par] << "  " << rbr_th_2b[par] << endl;
		rbr_th_sum+=rbr_th[par];
		rbr_th_2b_sum+=rbr_th_2b[par];
	}
	cout << setw(10) << "sum" << "  " << setw(10) << rbr_th_sum << "  " << rbr_th_2b_sum << endl;

	// for comparison 
	// nda relative br (2body)(use nominal)
	double rbr_nda_n = rbr_2b[1]/(rbr_2b[1]+rbr_2b[3]+rbr_2b[6]);
	double rbr_nda_da = (rbr_2b[3]+rbr_2b[6])/(rbr_2b[1]+rbr_2b[3]+rbr_2b[6]);
	cout << "Relative BR (n vs d/a): n   : " << setprecision(2) << rbr_nda_n << endl;
	cout << "Relative BR (n vs d/a): d/a : " << rbr_nda_da << endl;

	
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
		h_Ex_particle_th[p]->Scale(1.0/h_Ex->GetEntries());
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
	c_Ex->Print(pdfname.c_str());
	c_Ex->Update();
	c_Ex->Clear();
	//
	h_Ex->SetTitle("w/ energy th. (two-body decay)");
	h_Ex->GetXaxis()->SetRangeUser(8,47);
	h_Ex->Draw("HIST");
	for(int p=0;p<num_particle;p++){
		h_Ex_particle_th[p]->Draw("HISTsame");
	}
	if(Ex_max>8 && Ex_max<47) l_Ex_max->Draw("same");
	if(Ex_min>8 && Ex_min<47) l_Ex_min->Draw("same");
	c_Ex->Print(pdfname.c_str());
	c_Ex->Update();
	c_Ex->Clear();
	//
	THStack* h_s_Ex_th = new THStack("h_s_Ex_th","");
	h_Ex->Draw("HIST");
	for(int p=2;p<num_particle;p++){
		h_Ex_particle_th[p]->SetFillColor(color_root[p]);
		h_s_Ex_th->Add(h_Ex_particle_th[p]);
	}
	h_s_Ex_th->Draw("HISTsame");
	c_Ex->Print(pdfname.c_str());
	c_Ex->Update();
	c_Ex->Clear();
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
		if(kE_th[p]>0){
			TLine* l_kE_th = new TLine(kE_th[p],0,kE_th[p],h_kE[p]->GetMaximum()*1.05);
			l_kE_th->SetLineStyle(2);
			l_kE_th->Draw("same");
		}
		os.str("");
		os << "Prob = " << fixed << setprecision(2) << (double)h_kE[p]->GetEntries()/numofevent;
		TText* t_kE = new TText(5,h_kE[p]->GetMaximum()*0.7,os.str().c_str());
		t_kE->SetTextSize(0.06);
		t_kE->Draw("same");
		//
		os.str("");
		os << "Prob (w/ th)= " << fixed << (double)count_particle_th[p]/numofevent << " (" << (double)count_particle_th[p]/h_kE[p]->GetEntries()*100 << "%)";
		TText* t_kE_th = new TText(5,h_kE[p]->GetMaximum()*0.6,os.str().c_str());
		t_kE_th->SetTextSize(0.06);
		t_kE_th->Draw("same");
	}
	c_kE->cd(8);
	TPaveText* t = new TPaveText(0.1,0.1,0.9,0.9);
	os.str("");
	os << "# generated events = " << numofevent;
	t->AddText(os.str().c_str());
	for(int par=0;par<num_particle;par++){
		if(kE_th[par]>0){
			os.str("");
			os << "kE_th(" << particle_name[par].substr(0,1) << ") = " << setprecision(1) << kE_th[par];
			t->AddText(os.str().c_str());
		}
	}
	t->Draw("same");
	//
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_kE" << suffix.c_str() << ".pdf";
	c_kE->Print(os.str().c_str());


	//----------------------------//
	// nda plot
	const int nda_data=3;
	const int nda_color[nda_data]={416+2,800+6,600-6};
	const int	nda_color_this = 616-6;
	string nda_data_name[nda_data]={"SAbe_Tohoku_2023","Hu","Panin"};
	string nda_data_legend[nda_data]={"Abe et al. (TALYS)","Hu et al. (TALYS)","Exp. by Panin et al."};
	TFile* rootf_nda[nda_data];
	TH1D* h_nda[nda_data];
	for(int i=0;i<nda_data;i++){
		os.str("");
		os << "data/" << nda_data_name[i].c_str() << "/" << nda_data_name[i].c_str() << ".root";
		rootf_nda[i] = new TFile(os.str().c_str(),"READ");
		TEnv* env = (TEnv*) rootf_nda[i]->Get("env");
		//
		os.str("");
		os << "h_nda_" << nda_data_name[i].c_str();
		h_nda[i] = new TH1D(os.str().c_str(),"",20,-10,10);
		h_nda[i]->Fill(i,env->GetValue("rbr_nda_n",-9999.));
		h_nda[i]->Fill(i+nda_data+3,env->GetValue("rbr_nda_da",-9999.));
		h_nda[i]->SetFillColor(nda_color[i]);
		//h_nda[i]->SetLineWidth(2);
		if(i==2){
			 h_nda[i]->SetFillStyle(3445);
				h_nda[i]->SetLineColor(nda_color[i]);
		}
		else h_nda[i]->SetLineColor(1);
	}
	TH1D* h_nda_this = new TH1D("h_nda_this","",20,-10,10);
	h_nda_this->Fill(-1,rbr_nda_n*1e2);
	h_nda_this->Fill(-1+nda_data+3,rbr_nda_da*1e2);
	h_nda_this->SetFillColor(nda_color_this);
	h_nda_this->SetLineColor(1);


	TCanvas* c_nda =new TCanvas("c_nda","c_nda",0,0,800,600);
	gPad->SetRightMargin(0.02);
	gPad->SetTopMargin(0.02);
	TH1F* waku_nda = gPad->DrawFrame(-2,0,10,88);
	waku_nda->GetXaxis()->SetLabelSize(0);
	waku_nda->GetXaxis()->SetTickSize(0);
	waku_nda->GetYaxis()->SetTitle("Relative branching ratio (%)");
	waku_nda->GetYaxis()->CenterTitle();
	for(int i=0;i<nda_data;i++){
		h_nda[i]->Draw("HISTsame");
	}
	h_nda_this->Draw("HISTsame");
	TLegend* leg_nda = new TLegend(0.65,0.6,0.97,0.97);
	leg_nda->SetBorderSize(0);
	leg_nda->SetFillStyle(0);
	leg_nda->AddEntry(h_nda_this,"This work (TALYS)","F");
	for(int i=0;i<nda_data;i++){
		os.str("");
		os << nda_data_legend[i].c_str();
		leg_nda->AddEntry(h_nda[i],os.str().c_str(),"F");
	}
	leg_nda->Draw("same");
	TLatex* l_n = new TLatex(1,-8,"n");
	l_n->SetTextFont(12);
	l_n->SetTextSize(0.07);
	l_n->Draw("same");
	TLatex* l_da = new TLatex(0.5+nda_data+3,-8,"d/#alpha");
	l_da->SetTextFont(12);
	l_da->SetTextSize(0.07);
	l_da->Draw("same");
	gPad->RedrawAxis();
	//
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_nda" << suffix.c_str() << ".pdf";
	c_nda->Print(os.str().c_str());



	//---  br plot
	const int br_data=4;
	const int	br_color_this = 616-6;
	TH1D* h_br_this  = new TH1D("h_br_this","",50,-10,40);
	TH1D* h_br_2b_this  = new TH1D("h_br_2b_this","",50,-10,40);
	h_br_this->SetLineColor(br_color_this);
	//h_br_2b_this->SetFillStyle(3445);
	//h_br_2b_this->SetLineColor(br_color_this);
	h_br_2b_this->SetLineColor(1);
	h_br_2b_this->SetFillColor(br_color_this);
	//
	h_br_this->Fill(-1,rbr_th[1]/2);
	h_br_2b_this->Fill(-1,rbr_th_2b[1]/2);
	h_br_this->Fill(-1+br_data*1.5,rbr_th[2]);
	h_br_2b_this->Fill(-1+br_data*1.5,rbr_th_2b[2]);
	h_br_this->Fill(-1+br_data*3,rbr_th[3]);
	h_br_2b_this->Fill(-1+br_data*3,rbr_th_2b[3]);
	h_br_this->Fill(-1+br_data*4.5,rbr_th[4]);
	h_br_2b_this->Fill(-1+br_data*4.5,rbr_th_2b[4]);
	h_br_this->Fill(-1+br_data*6,rbr_th[6]);
	h_br_2b_this->Fill(-1+br_data*6,rbr_th_2b[6]);


	const int br_color[br_data]={416+2,800+6,600-7,1};
	string br_data_name[br_data]={"SAbe_Tohoku_2023","Hu","CASCADE","Yosoi"};
	string br_data_legend[br_data]={"Abe et al. (TALYS)","Hu et al. (TALYS)","Yosoi et al. (CASCADE)", "Exp. by Yosoi et al."};
	TFile* rootf_br[br_data];
	TH1D* h_br[br_data];
	TH1D* h_br_2b[br_data];
	const int bin_offset=11;
	double chi2 =0;
	for(int i=0;i<br_data;i++){
		os.str("");
		os << "data/" << br_data_name[i].c_str() << "/" << br_data_name[i].c_str() << ".root";
		if(br_data_name[i]=="CASCADE"){
			os.str("");
			os << "data/Yosoi/Yosoi.root";
		}
		cout << os.str().c_str() << endl;
		rootf_br[i] = new TFile(os.str().c_str(),"READ");
		TEnv* env = (TEnv*) rootf_br[i]->Get("env");
		//
		os.str("");
		os << "h_br_" << br_data_name[i].c_str();
		h_br[i] = new TH1D(os.str().c_str(),"",50,-10,40);
		os.str("");
		os << "h_br_2b_" << br_data_name[i].c_str();
		h_br_2b[i] = new TH1D(os.str().c_str(),"",50,-10,40);
		//
		if(br_data_name[i]=="CASCADE"){
			h_br[i]->SetBinContent(i+bin_offset,env->GetValue("cas_br_n",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset,env->GetValue("cas_br_n_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*1.5,env->GetValue("cas_br_p",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*1.5,env->GetValue("cas_br_p_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*3,env->GetValue("cas_br_d",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*3,env->GetValue("cas_br_d_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*4.5,env->GetValue("cas_br_t",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*4.5,env->GetValue("cas_br_t_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*6,env->GetValue("cas_br_a",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*6,env->GetValue("cas_br_a_2b",-9999.));
		}else{
			h_br[i]->SetBinContent(i+bin_offset,env->GetValue("br_n",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset,env->GetValue("br_n_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*1.5,env->GetValue("br_p",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*1.5,env->GetValue("br_p_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*3,env->GetValue("br_d",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*3,env->GetValue("br_d_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*4.5,env->GetValue("br_t",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*4.5,env->GetValue("br_t_2b",-9999.));
			h_br[i]->SetBinContent(i+bin_offset+br_data*6,env->GetValue("br_a",-9999.));
			h_br_2b[i]->SetBinContent(i+bin_offset+br_data*6,env->GetValue("br_a_2b",-9999.));
			if(br_data_name[i]=="Yosoi"){
				h_br[i]->SetBinError(i+bin_offset,env->GetValue("br_e_n",-9999.));
				h_br_2b[i]->SetBinError(i+bin_offset,env->GetValue("br_e_n_2b",-9999.));
				h_br[i]->SetBinError(i+bin_offset+br_data*1.5,env->GetValue("br_e_p",-9999.));
				h_br_2b[i]->SetBinError(i+bin_offset+br_data*1.5,env->GetValue("br_e_p_2b",-9999.));
				h_br[i]->SetBinError(i+bin_offset+br_data*3,env->GetValue("br_e_d",-9999.));
				h_br_2b[i]->SetBinError(i+bin_offset+br_data*3,env->GetValue("br_e_d_2b",-9999.));
				h_br[i]->SetBinError(i+bin_offset+br_data*4.5,env->GetValue("br_e_t",-9999.));
				h_br_2b[i]->SetBinError(i+bin_offset+br_data*4.5,env->GetValue("br_e_t_2b",-9999.));
				h_br[i]->SetBinError(i+bin_offset+br_data*6,env->GetValue("br_e_a",-9999.));
				h_br_2b[i]->SetBinError(i+bin_offset+br_data*6,env->GetValue("br_e_a_2b",-9999.));
				//
				chi2 += pow( (rbr_th[2] - env->GetValue("br_p",-9999.))/env->GetValue("br_e_p",-9999.), 2);
				chi2 += pow( (rbr_th_2b[2] - env->GetValue("br_p_2b",-9999.))/env->GetValue("br_e_p_2b",-9999.), 2);
				chi2 += pow( (rbr_th[3] - env->GetValue("br_d",-9999.))/env->GetValue("br_e_d",-9999.), 2);
				chi2 += pow( (rbr_th_2b[3] - env->GetValue("br_d_2b",-9999.))/env->GetValue("br_e_d_2b",-9999.), 2);
				chi2 += pow( (rbr_th[4] - env->GetValue("br_t",-9999.))/env->GetValue("br_e_t",-9999.), 2);
				chi2 += pow( (rbr_th_2b[4] - env->GetValue("br_t_2b",-9999.))/env->GetValue("br_e_t_2b",-9999.), 2);
				chi2 += pow( (rbr_th[6] - env->GetValue("br_a",-9999.))/env->GetValue("br_e_a",-9999.), 2);
				chi2 += pow( (rbr_th_2b[6] - env->GetValue("br_a_2b",-9999.))/env->GetValue("br_e_a_2b",-9999.), 2);
			}
		}
		//
		h_br[i]->SetLineColor(br_color[i]);
		//h_br_2b[i]->SetFillStyle(3445);
		if(i==3)  h_br_2b[i]->SetFillStyle(3445);
		//h_br_2b[i]->SetLineColor(br_color[i]);
		h_br_2b[i]->SetLineColor(1);
		h_br_2b[i]->SetFillColor(br_color[i]);
	}


	const double max_br_plot=34;
	TCanvas* c_br =new TCanvas("c_br","c_br",0,0,800,600);
	gPad->SetRightMargin(0.02);
	gPad->SetTopMargin(0.02);
	TH1F* waku_br = gPad->DrawFrame(-3,0,30,max_br_plot);
	waku_br->GetXaxis()->SetLabelSize(0);
	waku_br->GetXaxis()->SetTickSize(0);
	waku_br->GetYaxis()->SetTitle("Branching ratio (%)");
	waku_br->GetYaxis()->CenterTitle();
	h_br_this->Draw("HISTsame");
	h_br_2b_this->Draw("HISTsame");
	TLegend* leg_br = new TLegend(0.6,0.62,0.97,0.97);
	leg_br->SetBorderSize(0);
	leg_br->SetFillStyle(0);
	leg_br->AddEntry(h_br_2b_this,"This work (TALYS)","f");
	for(int i=0;i<br_data;i++){
		h_br[i]->Draw("HISTsame");
		h_br_2b[i]->Draw("HISTsame");
		if(br_data_name[i]=="Yosoi"){
			h_br[i]->Draw("Esame");
			h_br_2b[i]->Draw("Esame");
		}
		os.str("");
		os << br_data_legend[i].c_str();
		leg_br->AddEntry(h_br_2b[i],os.str().c_str(),"f");
	}
	leg_br->Draw("same");
	//
	TLine* line_br_n = new TLine(br_data,0,br_data,max_br_plot);
	line_br_n->SetLineStyle(2);
	line_br_n->Draw("same");
	TLatex* l_br_n_factor = new TLatex(0,max_br_plot*0.9,"#times 1/2");
	l_br_n_factor->SetTextSize(0.06);
	l_br_n_factor->Draw("same");
	os.str("");
	os << "#chi^{2} = " << chi2;
	TLatex* l_chi2 = new TLatex(br_data*1.5,max_br_plot*0.9,os.str().c_str());
	l_chi2->SetTextSize(0.07);
	//l_chi2->Draw("same");
	//
	int index=0;
	for(int p=1;p<num_particle;p++){
		if(particle_name[p]=="helium-3") continue;
		TLatex* l_br;
		if(particle_name[p]=="alpha") l_br = new TLatex(br_data*index*1.5+1,-2.5,"#alpha");
		else l_br = new TLatex(br_data*index*1.5+1,-2.5,particle_name[p].substr(0,1).c_str());
		l_br->Draw("same");
		l_br->SetTextFont(12);
		l_br->SetTextSize(0.07);
		index++;
	}
	gPad->RedrawAxis();
	//
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_br" << suffix.c_str() << ".pdf";
	c_br->Print(os.str().c_str());


	TCanvas* c_Ex_kE = new TCanvas("c_Ex_kE","c_Ex_kE",0,0,1200,600);
	c_Ex_kE->Divide(4,2);
	for(int p=0;p<num_particle;p++){
		c_Ex_kE->cd(p+1);
		h_Ex_kE[p]->SetTitle(particle_name[p].c_str());
		h_Ex_kE[p]->SetStats(0);
		h_Ex_kE[p]->GetXaxis()->SetRangeUser(0,50);
		h_Ex_kE[p]->GetYaxis()->SetRangeUser(0,28);
		h_Ex_kE[p]->GetXaxis()->SetTitle("Excitation energy (MeV)");
		h_Ex_kE[p]->GetYaxis()->SetTitle("Kinetic energy of particle (MeV)");
		h_Ex_kE[p]->Draw("colz");
		//TLine* l_Ex_kE = new TLine(SE[p],0,50,
		TF1* f_Ex_kE = new TF1("f_Ex_kE","(x-[1])*[0]",0,50);
		double factor=1;
		if(p==1 && p==2) factor= 1- 1.0/11;
		if(p==3) factor = 1-2.0/11;
		if(p==4 || p==5) factor = 1-3.0/11;
		if(p==6) factor = 1-4.0/11;
		f_Ex_kE->SetParameter(0,factor);
		f_Ex_kE->SetParameter(1,SE[p]);
		f_Ex_kE->Draw("same");
	}
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_Ex_kE" << suffix.c_str() << ".pdf";
	c_Ex_kE->Print(os.str().c_str());


	TCanvas* c_SE = new TCanvas("c_SE","c_SE",0,0,1200,600);
	c_SE->Divide(4,2);
	for(int p=0;p<num_particle;p++){
		c_SE->cd(p+1);
		h_SE[p]->SetTitle(particle_name[p].c_str());
		h_SE[p]->SetStats(0);
		h_SE[p]->GetXaxis()->SetRangeUser(0,50);
		h_SE[p]->GetXaxis()->SetTitle("Separation E + Excitation daughter (MeV)");
		h_SE[p]->Draw("HIST");
		TLine* l_SE = new TLine(SE[p],0,SE[p],h_SE[p]->GetMaximum()*1.05);
		l_SE->SetLineStyle(2);
		l_SE->Draw("same");
		//
		TLine* l_3b = new TLine(criteria_3b[p],0,criteria_3b[p],h_SE[p]->GetMaximum()*1.05);
		l_3b->SetLineStyle(2);
		l_3b->Draw("same");
		//
		os.str("");
		os << "Prob = " << (double)h_SE[p]->GetEntries()/numofevent*100;
		TText* t_SE  = new TText(20,h_SE[p]->GetMaximum()*0.85,os.str().c_str());
		t_SE->Draw("same");
		//
		os.str("");
		os << "Prob (3b) = " << (double)count_particle_3b[p]/numofevent*100;
		TText* t_SE_3b  = new TText(20,h_SE[p]->GetMaximum()*0.75,os.str().c_str());
		t_SE_3b->Draw("same");
	}
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_SE" << suffix.c_str() << ".pdf";
	c_SE->Print(os.str().c_str());


	// output 
	os.str("");
	os << "sim_out/Br_" << target.c_str() << "_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << ".root";
	TFile* outf = new TFile(os.str().c_str(),"RECREATE");
	TEnv env_o("env");
	env_o.SetValue("rbr_nda_n",rbr_nda_n*100);
	env_o.SetValue("rbr_nda_da",rbr_nda_da*100);
	env_o.SetValue("br_n_2b",rbr_th_2b[1]/2);
	env_o.SetValue("br_n",rbr_th[1]/2);
	env_o.SetValue("br_p_2b",rbr_th_2b[2]/2);
	env_o.SetValue("br_p",rbr_th[2]/2);
	env_o.SetValue("br_d_2b",rbr_th_2b[3]/2);
	env_o.SetValue("br_d",rbr_th[3]/2);
	env_o.SetValue("br_t_2b",rbr_th_2b[4]/2);
	env_o.SetValue("br_t",rbr_th[4]/2);
	env_o.SetValue("br_a_2b",rbr_th_2b[6]/2);
	env_o.SetValue("br_a",rbr_th[6]/2);
	env_o.Write("env");
	//
	for(int i=0;i<nda_data;i++){
		h_nda[i]->Write();
	}
	for(int i=0;i<br_data;i++){
		h_br[i]->Write();
		h_br_2b[i]->Write();
	}
	h_nda_this->Write();
	h_br_this->Write();
	h_br_2b_this->Write();
	//
	//outf->Close();
	//delete outf;

	return 0;
}
