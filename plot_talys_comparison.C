#include "include/consts.hh"
#include <sstream>

using namespace std;

int plot_talys_comparison(){
	// ---- FIXME ---- //
	string name_target="11C";
	//string name_target="11B";
  const float max_Ex_plot=50;
	// --------------- //


	ostringstream os;

	// New (TALYS 1.96 on this docker)
	TFile* rootf_new = new TFile( ((string)"output/Br_" + name_target + (string)".root").c_str(),"READ");
	TGraph* g_target_br[num_particle];
	for(int i=0;i<num_particle;i++){
		os.str("");
		os << "g_" << name_target.c_str() << "_br_" << i;
		g_target_br[i] = (TGraph*) rootf_new->Get(os.str().c_str());
	}
	

	// Old (TALYS 1.95 on komachi (tohoku rcns server))
	TFile* rootf_old = new TFile( ((string)"output_rcns/Br_" + name_target + (string)"_summary.root").c_str(),"READ");
	TGraph* g_old_target_br[num_particle];
	for(int i=0;i<num_particle;i++){
		os.str("");
		os << "g_" << name_target.c_str() << "_br_" << particle_name[i].substr(0,1);
		g_old_target_br[i] = (TGraph*) rootf_old->Get(os.str().c_str());
	}

	TLegend* leg_target_br = new TLegend(0.4,0.7,0.9,0.9);
	leg_target_br->SetNColumns(3);
	for(int p=0;p<num_particle;p++){
		leg_target_br->AddEntry(g_target_br[p],particle_name[p].c_str(),"PL");
	}
	
	gStyle->SetTitleXSize(0.045);
	gStyle->SetTitleYSize(0.045);
	gStyle->SetTitleYOffset(0.95);

	TCanvas* c_target_br = new TCanvas("c_target_br","",0,0,1600,600);
	c_target_br->Divide(2);
	//
	c_target_br->cd(1);
	TH1F* waku_old_target_br = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
	os.str("");
	os << "TALYS v1.95: " << name_target.c_str();
	waku_old_target_br->SetTitle(os.str().c_str());
	waku_old_target_br->GetXaxis()->SetTitle("Excitation energy [MeV]");
	waku_old_target_br->GetYaxis()->SetTitle("Branching ratio");
	gPad->SetGrid();
	for(int p=0;p<num_particle;p++){
		g_old_target_br[p]->Draw("PLsame");
	}
	leg_target_br->Draw("same");
	//
	c_target_br->cd(2);
	TH1F* waku_target_br = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
	os.str("");
	os << "TALYS v1.96: " << name_target.c_str();
	waku_target_br->SetTitle(os.str().c_str());
	waku_target_br->GetXaxis()->SetTitle("Excitation energy [MeV]");
	waku_target_br->GetYaxis()->SetTitle("Branching ratio");
	gPad->SetGrid();
	for(int p=0;p<num_particle;p++){
		g_target_br[p]->Draw("PLsame");
	}
	leg_target_br->Draw("same");
	gPad->RedrawAxis();
	//
	os.str("");
	os << "fig/" << "fig_" << name_target.c_str() << "_br_comparison.pdf";
	c_target_br->Print(os.str().c_str());




	TGraph* g_target_br_dev[num_particle];
	for(int p=0;p<num_particle;p++){
		g_target_br_dev[p] = new TGraph;
		os.str("");
		os << "g_target_br_dev_" << p;
		g_target_br_dev[p]->SetName(os.str().c_str());
		g_target_br_dev[p]->SetMarkerStyle(7);
		g_target_br_dev[p]->SetMarkerColor(color_root[p]);
		g_target_br_dev[p]->SetLineWidth(2);
		g_target_br_dev[p]->SetLineColor(color_root[p]);

		//
		int index=0;
		for(int k=0;k<g_target_br[p]->GetN();k++){
			double ex, br_new,br_old;
			g_target_br[p]->GetPoint(k,ex,br_new);
			br_old = g_old_target_br[p]->Eval(ex);

			double dev;
			if(br_old>0){
				dev = (br_new-br_old)/br_old*100;//(%)
			}else{
				continue;
			}
			g_target_br_dev[p]->SetPoint(index,ex,dev);
			index++;
		}
	}
	

	TCanvas* c_target_br_dev = new TCanvas("c_target_br_dev","",0,0,800,600);
	TH1F* waku_target_br_dev = gPad->DrawFrame(0,-100,max_Ex_plot,100);
	os.str("");
	os << "TALYS deviation based on v1.95: " << name_target.c_str();
	waku_target_br_dev->SetTitle(os.str().c_str());
	waku_target_br_dev->GetXaxis()->SetTitle("Excitation energy [MeV]");
	waku_target_br_dev->GetYaxis()->SetTitle("Deviation (%)");
	gPad->SetGrid();
	for(int p=0;p<num_particle;p++){
		g_target_br_dev[p]->Draw("PLsame");
	}
	//leg_target_br->Draw("same");
	os.str("");
	os << "fig/" << "fig_" << name_target.c_str() << "_br_dev.pdf";
	c_target_br_dev->Print(os.str().c_str());












	return 0;
}
