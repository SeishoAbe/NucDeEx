int plot_Ex(){
	//string target = "C";
	string target = "O";

	gStyle->SetTextFont(132);
  gStyle->SetTextSize(0.08);
  gStyle->SetTitleSize(0.06,"XYZ");
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetLabelSize(0.06,"XYZ");
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetLegendFont(132);
  //gStyle->SetLegendTextSize(0.04);
  gStyle->SetTitleYOffset(1.05);


	ostringstream os;
	os << "output_nuwro/histogram_deex_SF_LFG_14_1000_CCQE_" << target.c_str()
	   << ".root";
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	TH1D* h_Ex[4];
	for(int i=0;i<4;i++){
		os.str("");
		os << "h_Ex_" << i;
		h_Ex[i] = (TH1D*)rootf->Get(os.str().c_str());
		h_Ex[i]->Scale(1./h_Ex[0]->GetEntries());
	}


	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	gPad->SetRightMargin(0.05);
	gPad->SetTopMargin(0.05);
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.12);
	h_Ex[0]->GetXaxis()->CenterTitle();
	h_Ex[0]->GetYaxis()->CenterTitle();
	h_Ex[0]->GetXaxis()->SetTitle("Excitation energy (MeV)");
	h_Ex[0]->GetYaxis()->SetTitle("A.U");
	h_Ex[0]->GetXaxis()->SetTitleSize(0.06);
	h_Ex[0]->GetYaxis()->SetTitleSize(0.06);
	h_Ex[0]->GetYaxis()->SetNdivisions(505);
	h_Ex[0]->GetYaxis()->SetMaxDigits(4);
	for(int i=0;i<4;i++){
		h_Ex[i]->Draw("HISTsame");
	}
	os.str("");
	os << "#font[12]{p}_{1/2}-hole: " << fixed << setprecision(1) << (double)h_Ex[3]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
	TLatex* t_p12 = new TLatex(40,h_Ex[0]->GetMaximum()*0.5,os.str().c_str());
	if(target=="O") t_p12->Draw("same");
	os.str("");
	os << "#font[12]{p}_{3/2}-hole: " << fixed << setprecision(1) << (double)h_Ex[2]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
	TLatex* t_p32 = new TLatex(40,h_Ex[0]->GetMaximum()*0.35,os.str().c_str());
	t_p32->Draw("same");
	os.str("");
	os << "#font[12]{s}_{1/2}-hole: " << fixed << setprecision(1) << (double)h_Ex[1]->GetEntries()/h_Ex[0]->GetEntries()*100 << "%";
	TLatex* t_s12 = new TLatex(40,h_Ex[0]->GetMaximum()*0.2,os.str().c_str());
	t_s12->Draw("same");
	gPad->RedrawAxis();
	os.str("");
	os << "fig_nuwro/fig_Ex_paper_" << target.c_str() << ".pdf";
	c->Print(os.str().c_str());


	return 0;
}
