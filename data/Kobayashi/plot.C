
int plot(){
	TH1D* h = new TH1D("h","",26,3,8.20);
	TH1D* h_16_20 = new TH1D("h_16_20","",26,3,8.20);
	TH1D* h_20_30 = new TH1D("h_20_30","",26,3,8.20);
	TH1D* h_30_40 = new TH1D("h_30_40","",26,3,8.20);
	char buf[256];

	ifstream ifs_16_20("gamma_16_20.txt");
	while(ifs_16_20.getline(buf,sizeof(buf))){
		double e,ev;
		string st = buf;
		istringstream(st) >> e >> ev;
		h_16_20->Fill(e,ev);
		h->Fill(e,ev);
	}


	ifstream ifs_20_30("gamma_20_30.txt");
	while(ifs_20_30.getline(buf,sizeof(buf))){
		double e,ev;
		string st = buf;
		istringstream(st) >> e >> ev;
		h_20_30->Fill(e,ev);
		h->Fill(e,ev);
	}

	ifstream ifs_30_40("gamma_30_40.txt");
	while(ifs_30_40.getline(buf,sizeof(buf))){
		double e,ev;
		string st = buf;
		istringstream(st) >> e >> ev;
		h_30_40->Fill(e,ev);
		h->Fill(e,ev);
	}


	TCanvas* c_div = new TCanvas("c_div","c_div",0,0,1600,600);
	c_div->Divide(3);
	c_div->cd(1);
	h_16_20->SetMaximum(1700);
	h_16_20->SetMinimum(0);
	h_16_20->SetStats(0);
	h_16_20->Draw("HIST");
	c_div->cd(2);
	h_20_30->SetMaximum(1700);
	h_20_30->SetMinimum(0);
	h_20_30->SetStats(0);
	h_20_30->Draw("HIST");
	c_div->cd(3);
	h_30_40->SetMaximum(1700);
	h_30_40->SetMinimum(0);
	h_30_40->SetStats(0);
	h_30_40->Draw("HIST");
	c_div->Print("figure/fig_Eg.pdf");





	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	h->SetStats(0);
	h->Draw("HIST");
	TLine* l = new TLine(6,0,6,h->GetMaximum()*1.05);
	l->SetLineStyle(2);
	l->Draw("same");
	double ev_3_6=0, ev_6=0;
	for(int b=1;b<=26;b++){
		double e = h->GetBinCenter(b);
		//if(e>7.4) continue;
		if(e<=6) ev_3_6+=h->GetBinContent(b);
		else ev_6+=h->GetBinContent(b);
	}
	ostringstream os;
	os.str("");
	os << "Total: " << (int)(ev_3_6+ev_6) << " ev";
	TText* t_total = new TText(4.5,h->GetMaximum()*0.9,os.str().c_str());
	t_total->Draw("same");
	os.str("");
	os << "3 < E < 6 MeV: " << (int)ev_3_6 << " ev (" << fixed << setprecision(1) << ev_3_6/(ev_3_6+ev_6)*100 << "%)";
	TText* t_3_6 = new TText(4.5,h->GetMaximum()*0.8,os.str().c_str());
	t_3_6->Draw("same");
	os.str("");
	os << "E > 6 MeV: " << (int)ev_6 << " ev (" << fixed << setprecision(1) << ev_6/(ev_3_6+ev_6)*100 << "%)";
	TText* t_6 = new TText(4.5,h->GetMaximum()*0.7,os.str().c_str());
	t_6->Draw("same");
	c->Print("figure/fig_Eg_total.pdf");

	return 0;
}
