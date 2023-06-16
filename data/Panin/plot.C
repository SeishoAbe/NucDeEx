int plot(){
	gStyle->SetTextSize(0.08);
	gStyle->SetTitleSize(0.045);
	gStyle->SetTitleXSize(0.045);
	gStyle->SetTitleYSize(0.045);
	gStyle->SetTitleYOffset(0.95);

	char buf[256];

	TH1D* h_total = new TH1D("h_total","",50,0,50);
	double integ_total=0;
	ifstream ifs_total("scan_total.dat");
	while(ifs_total.getline(buf,sizeof(buf))){
		string st=buf;
		double x,y;
		istringstream(st) >> x >> y;
		if(y<0) y=0;
		h_total->Fill(x,y);
		if( x>16 && x<35) integ_total += y;
	}


	TH1D* h_n = new TH1D("h_n","",50,0,50);
	double integ_n=0;
	ifstream ifs_n("scan_n.dat");
	while(ifs_n.getline(buf,sizeof(buf))){
		string st=buf;
		double x,y;
		istringstream(st) >> x >> y;
		if(y<0) y=0;
		h_n->Fill(x,y);
		if( x>16 && x<35) integ_n += y;
	}

	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	h_total->SetStats(0);
	h_total->GetXaxis()->SetTitle("Excitatoin energy of ^{11}B (MeV)");
	h_total->GetYaxis()->SetTitle("A.U.");
	h_total->Draw("HIST");
	h_n->SetFillStyle(3445);
	h_n->SetFillColor(kBlue);
	h_n->Draw("HISTsame");
	TLine* l_min = new TLine(16,0,16,h_total->GetMaximum()*1.05);
	l_min->SetLineWidth(2);
	l_min->SetLineStyle(2);
	l_min->SetLineColor(kRed);
	l_min->Draw("same");
	TLine* l_max = new TLine(35,0,35,h_total->GetMaximum()*1.05);
	l_max->SetLineWidth(2);
	l_max->SetLineStyle(2);
	l_max->SetLineColor(kRed);
	l_max->Draw("same");
	c->Print("figure/fig_Ex.pdf");

	TFile* rootf = new TFile("Panin.root","RECREATE");
	h_n->Write();
	h_total->Write();
	TEnv env("env");
	double rbr_nda_n= integ_n/integ_total*100; // in (%)
	double rbr_nda_da= 100-rbr_nda_n;
	env.SetValue("rbr_nda_n",rbr_nda_n);
	env.SetValue("rbr_nda_da",rbr_nda_da);
	env.Print();
	env.Write("env");
	rootf->Close();
	delete rootf;


	return 0;
}
