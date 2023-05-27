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
	

	// prepare output root file
	os.str("");
	os << "output/Br_" << argv[1] << ".root";
	TFile* rootf = new TFile(os.str().c_str(),"RECREATE");


	// start Nucleus loop in NucleusTable
	for(int i=0;i<num_of_nuc;i++){
		if(!nucleus_table->GetNucleusPtr(i)->flag_data) continue; // no data -> continue

		Nucleus* nuc_target = nucleus_table->GetNucleusPtr(i); // get pointer of Nulceus
		string name_target = (string)nuc_target->name;

		// check info in Nucleus 
		if(! (nuc_target->CheckTotalPop() && nuc_target->CheckPop()) ){
			cerr << "ERROR: There was a bug in population info: " << nuc_target->name << endl;
			return 0;
		}else{
			cout << "Population for " << nuc_target->name << " looks OK" << endl;
		}
		if(! (nuc_target->CheckEx())) {
			cerr << "ERROR: There was unexpected behaviour in Ex" << endl;
			return 0;
		}
		int bin_target = nuc_target->Ex_bin[0];




		// Prepare TGraph
		// population TGraph
		TGraph* g_target_pop[parity];
		int index[parity]={0};
		for(int par=0;par<parity;par++){
			g_target_pop[par] = new TGraph;
			os.str("");
			os << "g_" << name_target.c_str() << "_pop_" << par;
			g_target_pop[par]->SetName(os.str().c_str());
		}
		// target br TGraph
		TGraph* g_target_br[num_particle];
		int	index_br[num_particle]={0};
		for(int p=0;p<num_particle;p++){
			index_br[p]=0;
			g_target_br[p] = new TGraph();
			os.str("");
			os << "g_" << name_target.c_str() << "_br_" << p;
			g_target_br[p]->SetName(os.str().c_str());
			g_target_br[p]->SetMarkerStyle(7);
			g_target_br[p]->SetMarkerColor(color_root[p]);
			g_target_br[p]->SetLineWidth(2);
			g_target_br[p]->SetLineColor(color_root[p]);
			if(p==0) g_target_br[p]->SetPoint(index_br[p],0,1);
			else     g_target_br[p]->SetPoint(index_br[p],0,0);
			index_br[p]++;
		}
		// target br & ex of daugther
		TGraph* g_target_br_ex[num_particle][bin_target];
		int index_br_ex[num_particle][bin_target]={0};
		for(int p=0;p<num_particle;p++){
			for(int i=0;i<bin_target;i++){ // target bin loop
				os.str("");
				os << "g_" << name_target << "_br_ex_" << p << "_" << i;
				g_target_br_ex[p][i] = new TGraph;
				g_target_br_ex[p][i]->SetName(os.str().c_str());
				g_target_br_ex[p][i]->SetMarkerStyle(7);
				g_target_br_ex[p][i]->SetMarkerColor(color_root[p]);
				g_target_br_ex[p][i]->SetLineWidth(2);
				g_target_br_ex[p][i]->SetLineColor(color_root[p]);
				index_br_ex[p][i]=0;
			}
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
				//if(Ex_target<nuc_target->min_S() ){
				if(Ex_target<nuc_target->min_S() 
						|| (!nuc_target->flag_pop_data[i] && i<=1)){ 
					if(p==0) g_target_br[p]->SetPoint(index_br[p],Ex_target,1);
					else g_target_br[p]->SetPoint(index_br[p],Ex_target,0);
					index_br[p]++;
				}else if(pop_target_sum>0 && nuc_target->flag_pop_data[i]){
					float population = nuc_target->GetPopDaughterBinSum(p,i);//  (particle, exbin)
					g_target_br[p]->SetPoint(index_br[p],Ex_target,population/pop_target_sum);
					index_br[p]++;
					for(int j=0;j<nuc_target->Ex_bin_p[p][i] && population>0;j++){ // daughter ex bin loop
						float br_ex = nuc_target->pop_p[p][i][j]/population;
						//if(br_ex==0) cout << i << " " << p << " " << j << " " << br_ex << endl;
						g_target_br_ex[p][i]->SetPoint(index_br_ex[p][i],
																					 nuc_target->Ex_p[p][i][j],br_ex);
						index_br_ex[p][i]++;
					}
				}
			}
		}


		// --- Draw --- // 
		gStyle->SetTextSize(0.08);
		gStyle->SetTitleSize(0.045);
		gStyle->SetTitleXSize(0.045);
		gStyle->SetTitleYSize(0.045);
		gStyle->SetTitleYOffset(0.95);
		
		// --- Canvas --- //
		TCanvas* c_target_pop = new TCanvas("c_target_pop","",0,0,1200,600);
		c_target_pop->Divide(2);
		for(int par=0;par<parity;par++){
			c_target_pop->cd(par+1);
			TH1F* waku_target_pop = gPad->DrawFrame(0,0,max_Ex_plot,max_total_pop[par]*1.1);
			os.str("");
			os << "Population (normarized), " << name_target.c_str();
			waku_target_pop->SetTitle(os.str().c_str());
			waku_target_pop->GetXaxis()->SetTitle("Excitation energy [MeV]");
			waku_target_pop->GetYaxis()->SetTitle("Population (A.U.)");
			gPad->SetGrid();
			g_target_pop[par]->SetLineWidth(2);
			g_target_pop[par]->SetMarkerStyle(7);
			g_target_pop[par]->Draw("PLsame");
		}
		gPad->RedrawAxis();
		os.str("");
		os << "fig/" << argv[1] << "/fig_" << name_target.c_str() << "_pop.pdf";
		c_target_pop->Print(os.str().c_str());
		c_target_pop->Update();
		c_target_pop->Clear();
		delete c_target_pop;



		TCanvas* c_target_br = new TCanvas("c_target_br","",0,0,800,600);
		TH1F* waku_target_br = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
		os.str("");
		os << "Branching ratio: target: " << name_target.c_str();
		waku_target_br->SetTitle(os.str().c_str());
		waku_target_br->GetXaxis()->SetTitle("Excitation energy [MeV]");
		waku_target_br->GetYaxis()->SetTitle("Branching ratio");
		gPad->SetGrid();
		for(int p=0;p<num_particle;p++){
			g_target_br[p]->Draw("PLsame");
		}
		TLegend* leg_target_br = new TLegend(0.4,0.7,0.9,0.9);
		leg_target_br->SetNColumns(3);
		for(int p=0;p<num_particle;p++){
			leg_target_br->AddEntry(g_target_br[p],particle_name[p].c_str(),"PL");
		}
		leg_target_br->Draw("same");
		gPad->RedrawAxis();
		//
		os.str("");
		os << "fig/" << argv[1] << "/fig_" << name_target.c_str() << "_br.pdf";
		c_target_br->Print(os.str().c_str());
		c_target_br->Update();
		c_target_br->Clear();
		delete c_target_br;


		TCanvas* c_target_br_ex = new TCanvas("c_target_br_ex","",0,0,1200,600);
		string pdfname= (string)"fig/" + (string)argv[1] + (string)"/fig_" + name_target + (string)"_target_br_ex.pdf";
		c_target_br_ex->Print( (pdfname+(string)"[").c_str() );
		c_target_br_ex->Update();
		c_target_br_ex->Clear();
		for(int i=0;i<bin_target;i++){ // target bin loop
			if(nuc_target->Ex[0][i]<nuc_target->min_S()) continue;
			c_target_br_ex->Divide(4,2);
			for(int p=0;p<num_particle;p++){ // particle loop
				c_target_br_ex->cd(p+1);
				gPad->SetGrid();
				TH1F* waku = gPad->DrawFrame(0,0,max_Ex_plot,1.05);
				os.str("");
				os << particle_name[p];
				waku->SetTitle(os.str().c_str());
				waku->GetXaxis()->SetTitle("Excitation energy of daughter nucleus");
				waku->GetYaxis()->SetTitle("Relative branching ratio");
				if(g_target_br_ex[p][i]->GetN()==0){
					TText* text = new TText(max_Ex_plot*0.25,0.8,"NO DATA");
					text->SetTextColor(kGray+1);
					text->Draw("same");
				}else{
					g_target_br_ex[p][i]->Draw("PLsame");
					os.str("");
					os << "BR: " << scientific << setprecision(2) 
						 << nuc_target->GetPopDaughterBinSum(p,i)/nuc_target->GetPopParitySum(i);
					TText* text = new TText(max_Ex_plot*0.25,0.8,os.str().c_str());
					text->Draw("same");
				}
			}
			c_target_br_ex->cd(8);
			os.str("");
			os << name_target.c_str() << ", Ex[" << i << "] = " << fixed << nuc_target->Ex[0][i] << " (MeV)";
			TText* text_Ex = new TText(0.1,0.9,os.str().c_str());
			text_Ex->Draw("same");
			//
			gPad->RedrawAxis();
			os.str("");
			if(i==bin_target-1) os << (pdfname+(string)"]").c_str();
			else os << pdfname.c_str();
			c_target_br_ex->Print(os.str().c_str());
			c_target_br_ex->Update();
			c_target_br_ex->Clear();
		}
		delete c_target_br_ex;


		//
		rootf->cd();
		for(int par=0;par<parity;par++){
			g_target_pop[par]->Write();
		}
		for(int p=0;p<num_particle;p++){
			g_target_br[p]->Write();
			for(int i=0;i<bin_target;i++){ // target bin loop
				g_target_br_ex[p][i]->Write();
			}
		}
	}


	rootf->Close();
	delete rootf;


	delete read_talys;
	delete nucleus_table;


  return 0;
}
