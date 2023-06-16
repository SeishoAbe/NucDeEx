#include "include/consts.hh"
#include <TStyle.h>
#include <map>

using namespace std;

const double kE_th[num_particle]
	= {0,3.1,3.1,
		 4.0,4.6,0,4.5};

int plot_simulation_comparison(){
	// ---- FIXME ---- //
	string target = "11B";
	const double Ex_min =16;
	const double Ex_max =35;
		// negative -> not applied
	const bool flag_decay=1;
		// 0 -> use string w/ "g" (gamma)
		// 1- > use string w/o "g" (gamma) <- use this 
	// --------------- //

	ostringstream os;
	os.str("");
	os << "sim_out/" << target.c_str() << ".root";
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
	double rbr[num_particle]={0}; // two-body
	double rbr_2b[num_particle]={0}; // two-body
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);

		// --- Ex selection --- //
		h_Ex->Fill(Ex);
		if(Ex_min>0 && Ex<Ex_min) continue;
		if(Ex_max>0 && Ex>Ex_max) continue;

		if(max_size<size) max_size=size;
		int nmulti=0;
		int particle_counter=0; // w/o g
		int first_p=-1;
		for(int b=0;b<size;b++){
			if(PDG[b]==2112) nmulti++;
			int p=-1;
			for(int par=0;par<num_particle;par++){
				if(PDG_particle[par]==PDG[b]){
					p=par;
					break;
				}
			}
			if(p>=0){ // this is particle info (not daughter nuclei)
				h_kE[p]->Fill(kE[b]);
				if(p!=0) particle_counter++;
				if(b==0){
					first_p = p;
					rbr[p]++;
				}
			}
		}
		if(particle_counter==1) rbr_2b[first_p]++;
		itr = br.find(decay->c_str());
		if(itr != end(br) ) {
			itr->second += 1;
		}else{
			//cout << decay->c_str() << endl;
			br.insert(make_pair(decay->c_str(),1));
		}
		h_nmulti->Fill(nmulti);
		numofevent++;
	}
	cout << "max_size=" << max_size << endl;
	cout << "numofevent=" << numofevent << endl;

	double rbr_sum=0, rbr_2b_sum=0;
	for(int par=0;par<num_particle;par++){
		rbr[par] = rbr[par]/numofevent*100; // %
		rbr_2b[par] = rbr_2b[par]/numofevent*100; // %
		cout << setw(10) << particle_name[par].c_str() << "  " 
				 << setw(10) << rbr[par] << "  " << rbr_2b[par] << endl;
		rbr_sum+=rbr[par];
		rbr_2b_sum+=rbr_2b[par];
	}
	cout << setw(10) << "sum" << "  " << setw(10) << rbr_sum << "  " << rbr_2b_sum << endl;

	// for comparison 
	// nda relative br (2body)
	double rbr_nda_n = rbr_2b[1]/(rbr_2b[1]+rbr_2b[3]+rbr_2b[6]);
	double rbr_nda_da = (rbr_2b[3]+rbr_2b[6])/(rbr_2b[1]+rbr_2b[3]+rbr_2b[6]);
	cout << "Relative BR (n vs d/a): n   : " << rbr_nda_n << endl;
	cout << "Relative BR (n vs d/a): d/a : " << rbr_nda_da << endl;

	
	// output  (txt)
	os.str("");
	if(Ex_min>0) os << "_Exmin" << Ex_min;
	if(Ex_max>0) os	<< "_Exmax" << Ex_max;
	string suffix = os.str();
	os.str("");
	os << "sim_out/" << target.c_str() << suffix.c_str() << "_raw.txt";
	ofstream ofs_raw (os.str().c_str());
	os.str("");
	os << "sim_out/" << target.c_str() << suffix.c_str() << ".txt";
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
	}

	TCanvas* c_Ex = new TCanvas("c_Ex","",0,0,800,600);
	h_Ex->GetXaxis()->SetRangeUser(8,47);
	h_Ex->SetStats(0);
	h_Ex->SetMinimum(0);
	h_Ex->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_Ex->GetYaxis()->SetTitle("A.U.");
	h_Ex->Draw("HIST");
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
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_Ex" << suffix.c_str() << ".pdf";
	c_Ex->Print(os.str().c_str());


	TCanvas* c_Ex_1 = new TCanvas("c_Ex_1","",0,0,800,600);
	h_Ex->GetXaxis()->SetRangeUser(0,100);
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
	os.str("");
  //TF1* f_lorentzian = new TF1("f_lorentziaon","[0]/(TMath::Pi()*(TMath::Power(x-[1],2)+TMath::Power([0],2)))",0,100);
  //f_lorentzian->SetParameter(0,7); // HWHM
  //f_lorentzian->SetParameter(1,23); // mean
	//h_Ex->Fit(f_lorentzian,"","",18,30);
	//f_lorentzian->Draw("same");
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_Ex_1" << suffix.c_str() << ".pdf";
	c_Ex_1->Print(os.str().c_str());

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
	os << "fig_sim/fig_" << target.c_str() << "_nmulti" << suffix.c_str() << ".pdf";
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
	os << "fig_sim/fig_" << target.c_str() << "_kE" << suffix.c_str() << ".pdf";
	c_kE->Print(os.str().c_str());


	//----------------------------//
	// nda plot
	const int nda_data=3;
	const int nda_color[nda_data]={416+2,800+6,600-6};
	const int	nda_color_this = 616-6;
	string nda_data_name[nda_data]={"KamLAND","Hu","Panin"};
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
	}
	TH1D* h_nda_this = new TH1D("h_nda_this","",20,-10,10);
	h_nda_this->Fill(-1,rbr_nda_n*1e2);
	h_nda_this->Fill(-1+nda_data+3,rbr_nda_da*1e2);


	TCanvas* c_nda =new TCanvas("c_nda","c_nda",0,0,800,600);
	TH1F* waku_nda = gPad->DrawFrame(-2,0,10,85);
	waku_nda->GetXaxis()->SetLabelSize(0);
	waku_nda->GetXaxis()->SetTickSize(0);
	waku_nda->GetYaxis()->SetTitle("Relative branching ratio (%)");
	for(int i=0;i<nda_data;i++){
		h_nda[i]->Draw("HISTsame");
	}
	h_nda_this->SetFillColor(nda_color_this);
	h_nda_this->Draw("HISTsame");
	TLegend* leg = new TLegend(0.6,0.6,0.9,0.89);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(h_nda_this,"This work","F");
	for(int i=0;i<nda_data;i++){
		os.str("");
		os << nda_data_name[i].c_str();
		if(i!=0) os << " et al.";
		leg->AddEntry(h_nda[i],os.str().c_str(),"F");
	}
	leg->Draw("same");
	TLatex* l_n = new TLatex(1,-8,"n");
	l_n->SetTextFont(12);
	l_n->Draw("same");
	TLatex* l_da = new TLatex(0.5+nda_data+3,-8,"d/#alpha");
	l_da->SetTextFont(12);
	l_da->Draw("same");
	//
	os.str("");
	os << "fig_sim/fig_" << target.c_str() << "_nda" << suffix.c_str() << ".pdf";
	c_nda->Print(os.str().c_str());




	return 0;
}
