#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "consts.hh"
#include "NucleusTable.hh"
#include "ReadTALYS.hh"

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>

using namespace std;

int main(int argc, char* argv[]){
	if(argc!=2){
		cerr << argv[0] << " [Target nucleus]" << endl;
		return 0;
	}
  // --- FIXME  --- //
  const int Ex_max=60000; // keV
  const int Ex_bin_width=100; // keV
  const float max_Ex_plot=50;
  // ---------------//
  
  std::ostringstream os;
  NucleusTable* nucleus_table = new NucleusTable();
  if(!nucleus_table->ReadTables()){
		cerr << "something wrong" << endl;
		return 1;
	}
	os.str("");
	os << "output/output_" << argv[1];
	ReadTALYS* read_talys = new ReadTALYS(os.str().c_str(), nucleus_table);
	if(!read_talys->Read()){
		cerr << "something wrong happend..." << endl;
		return 0;
	}

	// Get general information
	//
	const int num_of_nuc = nucleus_table->GetNumofNuc();
	int num_nuc_data=0;
	for(int i=0;i<num_of_nuc;i++){
		if(nucleus_table->GetNucleusPtr(i)->flag_data==1)  num_nuc_data++;
	}
	cout << "# nuclei that has population data = " << num_nuc_data << endl;
	Nucleus* nuc_target = nucleus_table->GetNucleusPtr(argv[1]);
	// check info in Nucleus 
	if(! (nuc_target->CheckTotalPop() && nuc_target->CheckPop()) ){
		cerr << "ERROR: There was a bug in population info: " << nuc_target->name << endl;
		return 0;
	}
	if(! (nuc_target->CheckEx())) {
		cerr << "ERROR: There was unexpected behaviour in Ex" << endl;
		return 0;
	}
	int bin_target = nuc_target->Ex_bin[0];


	// ---- Draw --- // 
	gStyle->SetTitleSize(0.045);
	//
	TGraph* g_target_pop[parity];
	int index[parity]={0};
	for(int par=0;par<parity;par++){
		g_target_pop[par] = new TGraph;
		os.str("");
		os << "g_target_pop_" << par;
		g_target_pop[par]->SetName(os.str().c_str());
	}
	TGraph* g_target_br[num_particle];
	int	index_br[num_particle]={0};
	TGraph* g_target_br_ex[num_particle][bin_target];


	for(int p=0;p<num_particle;p++){
		g_target_br[p] = new TGraph();// ((string)"g_target_br_" + particle_name[p]).c_str());
		g_target_br[p]->SetMarkerStyle(7);
		//g_target_br[p]->SetMarkerSize(2);
		g_target_br[p]->SetMarkerColor(color_root[p]);
		g_target_br[p]->SetLineWidth(2);
		g_target_br[p]->SetLineColor(color_root[p]);
		if(p==0) g_target_br[p]->SetPoint(index_br[p],0,1);
		else     g_target_br[p]->SetPoint(index_br[p],0,0);
		index_br[p]++;
	}

	float max_total_pop[parity]={0};
	for(int i=0;i<bin_target;i++){ // target bin loop
		float Ex_target = nuc_target->Ex[0][i];
		float pop_target_sum = nuc_target->GetPopParitySum(i);

		for(int par=0;par<parity;par++){
			float pop_target_norm = nuc_target->pop[par][i]/nuc_target->total_pop[par];
			// fill normalized population
			g_target_pop[par]->SetPoint(index[par],Ex_target,pop_target_norm);
			if(max_total_pop[par]<pop_target_norm) max_total_pop[par] = pop_target_norm;
			index[par]++;
			
		}
		
		// calculate br
		for(int p=0;p<num_particle;p++){ // particle loop
			if(Ex_target<nuc_target->min_S()){
				if(p==0) g_target_br[p]->SetPoint(index_br[p],Ex_target,1);
				else g_target_br[p]->SetPoint(index_br[p],Ex_target,0);
				index_br[p]++;
			}else if(pop_target_sum>0){
				float population = nuc_target->GetPopDaughterBinSum(p,i);//  (particle, exbin)
				g_target_br[p]->SetPoint(index_br[p],Ex_target,population/pop_target_sum);
				index_br[p]++;
			}
		}
	}



	TCanvas* c_target_pop = new TCanvas("c_target_pop","",0,0,1200,600);
	c_target_pop->Divide(2);
	for(int par=0;par<parity;par++){
		c_target_pop->cd(par+1);
		TH1F* waku_target_pop = c_target_pop->DrawFrame(0,0,max_Ex_plot,max_total_pop[par]*1.1);
		os.str("");
		os << "Population (normarized), " << argv[1];
		waku_target_pop->SetTitle(os.str().c_str());
		waku_target_pop->GetXaxis()->SetTitle("Excitation energy [MeV]");
		waku_target_pop->GetYaxis()->SetTitle("Population (A.U.)");
		gPad->SetGrid();
		g_target_pop[par]->SetLineWidth(2);
		g_target_pop[par]->SetMarkerStyle(7);
		g_target_pop[par]->Draw("PLsame");
	}
	os.str("");
	os << "fig/fig_target_pop.pdf";
	c_target_pop->Print(os.str().c_str());
	c_target_pop->Update();
	c_target_pop->Clear();
	delete c_target_pop;

	TCanvas* c_target_br = new TCanvas("c_target_br","",0,0,800,600);
	TH1F* waku_target_br = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
	os.str("");
	os << "Branching ratio: target: " << argv[1];
	waku_target_br->SetTitle(os.str().c_str());
	waku_target_br->GetXaxis()->SetTitle("Excitation energy [MeV]");
	waku_target_br->GetYaxis()->SetTitle("Branching ratio");
	gPad->SetGrid();
	for(int p=0;p<num_particle;p++){
		g_target_br[p]->Draw("PLsame");
	}
	os.str("");
	os << "fig/fig_target_br.pdf";
	c_target_br->Print(os.str().c_str());
	c_target_br->Update();
	c_target_br->Clear();
	delete c_target_br;


	//
	os.str("");
	os << "output/Br_" << argv[1] << ".root";
	TFile* rootf = new TFile(os.str().c_str(),"RECREATE");
	for(int par=0;par<parity;par++){
		g_target_pop[par]->Write();
	}
	rootf->Close();
	delete rootf;


	delete read_talys;
	delete nucleus_table;


  return 0;
}
