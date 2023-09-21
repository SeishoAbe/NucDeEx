using namespace std;


int plot_comparison(){
	string target = "C";
	//string target = "O";
	
  gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.06);
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetLabelSize(0.06,"XYZ");
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetLegendFont(132);
  gStyle->SetLegendTextSize(0.045);
  gStyle->SetTitleYOffset(0.95);

	const int numfile=4; // FIXME//
	const int MDLQE=422;
	const double maxy=0.76;
	const double maxy_mean=2.3;
	TFile* rootf[numfile][2];
	// 1st: neut (this work)
	//		  neut (original deex)
	//      nuwro
	//      genie
	// 2nd: numu / numub
	TH1D* h_nmulti[numfile][2];
	ostringstream os;

	const string flavor[2]={"numu","numub"};
	const int pdg[2]={14,-14};

	for(int i=0;i<2;i++){

		// neut
		os.str("");	
		os << "output_neut/histogram_deex_neut_1GeV_" << flavor[i].c_str() << "_CCQE_" 
			 << target.c_str() << "_MDLQE" << MDLQE;
		if(target=="O")os<< "_NUCDEXITE0.root";
		else os <<".root";
		rootf[0][i] = new TFile(os.str().c_str(),"READ");
		cout << os.str().c_str() << endl;

		if(target=="O"){
			os.str("");	
			os << "output_neut/histogram_neut_1GeV_" << flavor[i].c_str() << "_CCQE_" 
				<< target.c_str() << "_MDLQE" << MDLQE << ".root";
			rootf[1][i] = new TFile(os.str().c_str(),"READ");
		  cout << os.str().c_str() << endl;
		}
		
		// nuwro
		os.str("");	
		os << "output_nuwro/histogram_deex_SF_LFG_" << pdg[i] << "_1000_CCQE_" 
			 << target.c_str() << ".root";
		rootf[2][i] = new TFile(os.str().c_str(),"READ");
		cout << os.str().c_str() << endl;

		// genie
		os.str("");	
		os << "output_genie/histogram_deex_CCQE." << pdg[i] << "." << target.c_str() 
		   << ".G18_10b_02_11a.root";
		rootf[3][i] = new TFile(os.str().c_str(),"READ");
		cout << os.str().c_str() << endl;
	}



	const int color[numfile]={416+1,416+1,600,616};
	const int style[numfile]={1,9,7,5};
	for(int f=0;f<numfile;f++){
		if(target!="O" && f==1) continue;
		for(int i=0;i<2;i++){
			if(f==1) 
				h_nmulti[f][i] = (TH1D*) rootf[f][i]->Get("h_nmulti_postFSI");
			else 
				h_nmulti[f][i] = (TH1D*) rootf[f][i]->Get("h_nmulti_postdeex");
			h_nmulti[f][i]->Scale(1./h_nmulti[f][i]->Integral());
			//h_nmulti[f][i]->Scale(1./h_nmulti[f][i]->GetEntries());
			h_nmulti[f][i]->SetLineColor(color[f]);
			h_nmulti[f][i]->SetLineStyle(style[f]);
			h_nmulti[f][i]->SetLineWidth(2);
		}
	}

		
	
	const double width=0.02;

	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	TH1F* waku = gPad->DrawFrame(-1,0,1,maxy_mean);
	waku->GetXaxis()->SetLabelSize(0);
	waku->GetXaxis()->SetTickLength(0);
	waku->GetXaxis()->CenterTitle();
	waku->GetYaxis()->CenterTitle();
	waku->GetYaxis()->SetTitle("Mean neutron multiplicity");
	TLine* l[numfile][2];
	TLegend* leg1 = new TLegend(0.55,0.15,0.95,0.5);
	for(int f=0;f<numfile;f++){
		if(target!="O" && f==1) continue;
		for(int i=0;i<2;i++){
			if(i==0) l[f][i] = new TLine(-0.8,h_nmulti[f][i]->GetMean(),-0.2,h_nmulti[f][i]->GetMean());
			else l[f][i] = new TLine(0.2,h_nmulti[f][i]->GetMean(),0.8,h_nmulti[f][i]->GetMean());
			l[f][i]->SetLineColor(color[f]);
			l[f][i]->SetLineStyle(style[f]);
			l[f][i]->SetLineWidth(2);
			l[f][i]->Draw("same");
			cout << h_nmulti[f][i]->GetMean() << endl;
		}
	}
	TLatex* l_nu = new TLatex(-0.55,-0.22,"#nu_{#mu}");
	l_nu->SetTextSize(0.09);
	l_nu->Draw("same");
	TLatex* l_nub = new TLatex(0.45,-0.22,"#bar{#nu}_{#mu}");
	l_nub->SetTextSize(0.09);
	l_nub->Draw("same");
	//
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0);
	leg1->AddEntry(h_nmulti[0][0],"NEUT (5.6.3)","l");
	if(target=="O") leg1->AddEntry(h_nmulti[1][0],"NEUT original deex (5.6.3)","l");
	leg1->AddEntry(h_nmulti[2][0],"NuWro (21.09.02)","l");
	leg1->AddEntry(h_nmulti[3][0],"GENIE (3.04.00)","l");
	leg1->Draw("same");
	os.str("");
	os << "fig_comparison/fig_deex_mean_nmulti_" << target.c_str() << ".pdf";
	c->Print(os.str().c_str());

	TCanvas* c_nu = new TCanvas("c_nu","c_nu",0,0,800,600);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	TLegend* leg = new TLegend(0.5,0.5,0.95,0.95);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry(h_nmulti[0][0],"NEUT (5.6.3)","l");
	if(target=="O") leg->AddEntry(h_nmulti[1][0],"NEUT original deex (5.6.3)","l");
	leg->AddEntry(h_nmulti[2][0],"NuWro (21.09.02)","l");
	leg->AddEntry(h_nmulti[3][0],"GENIE (3.04.00)","l");

	for(int f=0;f<numfile;f++){
		if(target!="O" && f==1) continue;
		h_nmulti[f][0]->SetStats(0);
		h_nmulti[f][0]->GetXaxis()->SetNdivisions(10);
		h_nmulti[f][0]->GetXaxis()->SetRangeUser(-0.5,6.5);
		h_nmulti[f][0]->GetYaxis()->SetTitle("A.U.");
		h_nmulti[f][0]->GetXaxis()->SetTitle("Neutron multiplicity");
		h_nmulti[f][0]->GetYaxis()->SetMaxDigits(3);
		h_nmulti[f][0]->SetMaximum(maxy);
		h_nmulti[f][0]->GetXaxis()->CenterTitle(1);
		h_nmulti[f][0]->GetYaxis()->CenterTitle(1);
		h_nmulti[f][0]->GetXaxis()->SetTitleSize(0.06);
		h_nmulti[f][0]->GetYaxis()->SetTitleSize(0.06);
		h_nmulti[f][0]->Draw("HISTsame");
	}
	leg->Draw("same");
	os.str("");
	os << "fig_comparison/fig_deex_nmulti_" << target.c_str() << "_nu.pdf";
	c_nu->Print(os.str().c_str());



	TCanvas* c_nub = new TCanvas("c_nub","c_nub",0,0,800,600);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	for(int f=0;f<numfile;f++){
		if(target!="O" && f==1) continue;
		h_nmulti[f][1]->SetStats(0);
		h_nmulti[f][1]->GetXaxis()->SetNdivisions(10);
		h_nmulti[f][1]->GetXaxis()->SetRangeUser(-0.5,6.5);
		h_nmulti[f][1]->GetYaxis()->SetTitle("A.U.");
		h_nmulti[f][1]->GetXaxis()->SetTitle("Neutron multiplicity");
		h_nmulti[f][1]->GetYaxis()->SetMaxDigits(3);
		h_nmulti[f][1]->SetMaximum(maxy);
		h_nmulti[f][1]->GetXaxis()->CenterTitle(1);
		h_nmulti[f][1]->GetYaxis()->CenterTitle(1);
		h_nmulti[f][1]->GetXaxis()->SetTitleSize(0.06);
		h_nmulti[f][1]->GetYaxis()->SetTitleSize(0.06);
		h_nmulti[f][1]->Draw("HISTsame");
	}
	leg->Draw("same");
	os.str("");
	os << "fig_comparison/fig_deex_nmulti_" << target.c_str() << "_nub.pdf";
	c_nub->Print(os.str().c_str());


	return 0;
}
