int plot_br_paper(){
	string nucleus="11B";
	string nuc_name="^{11}B^{*}";
	//string nucleus="11C";
	//string nuc_name="^{11}C^{*}";
	//string nucleus="15O";
	//string nuc_name="^{15}O^{*}";
	//string nucleus="15N";
	//string nuc_name="^{15}N^{*}";
	string path;
	if(nucleus.find("11")!=string::npos){
		path = "12C";
	}else{
		path = "16O";
	}
	const int bin=22;

	gStyle->SetTextFont(132);
	gStyle->SetTextSize(0.08);
	gStyle->SetTitleSize(0.06,"XYZ");
	gStyle->SetTitleFont(132,"XYZ");
	gStyle->SetLabelSize(0.06,"XYZ");
	gStyle->SetLabelFont(132,"XYZ");
	gStyle->SetLegendFont(132);
	gStyle->SetLegendTextSize(0.06);
	gStyle->SetTitleYOffset(0.90);

	ostringstream os;
	os.str("");
	os << "output/" << path.c_str() << "/Br_" << nucleus.c_str() << "_ldmodel2_parity_optmodall.root";
	TFile* rootf = new TFile(os.str().c_str(),"READ");
	const int part_num=7;
	const char* part_name[]={"#gamma","n","p","d","t","#font[132]{^{3}He}","#alpha"};
	TGraph* g[part_num];
	TGraph* g_ex[part_num];
	for(int i=0;i<part_num;i++){
		os.str("");
		os << "g_" << nucleus.c_str() << "_br_" << i;
		g[i] = (TGraph*) rootf->Get(os.str().c_str());
		//
		os.str("");
		os << "g_" << nucleus.c_str() << "_br_ex_" << i << "_" << bin;
		g_ex[i] = (TGraph*) rootf->Get(os.str().c_str());
	}
	const int order[part_num]={0,1,2,5,6,3,4};

	TCanvas* c = new TCanvas("c","c",0,0,800,600);
	gPad->SetRightMargin(0.02);
	gPad->SetTopMargin(0.02);
	gPad->SetLeftMargin(0.13);
	gPad->SetBottomMargin(0.13);
	TH1F* waku = gPad->DrawFrame(0,0,50,1.05);
	waku->GetXaxis()->SetTitle("Excitation energy (MeV)");
	waku->GetXaxis()->CenterTitle();
	waku->GetYaxis()->SetTitle("Branching ratio");
	waku->GetYaxis()->CenterTitle();
	TLegend* leg = new TLegend(0.6,0.65,0.97,0.97);
	leg->SetBorderSize(0);
	leg->SetTextFont(12);
	leg->SetFillStyle(0);
	leg->SetNColumns(4);
	for(int i=0;i<part_num;i++){
		g[i]->Draw("L");
	}
	for(int i=0;i<part_num;i++){
		leg->AddEntry(g[order[i]],part_name[order[i]],"L");
	}
	leg->Draw("same");
	TLatex* l_nucleus = new TLatex(18,0.85,nuc_name.c_str());
	l_nucleus->Draw("same");
	//l_nucleus->SetTextSize(0.06);
	os.str("");
	os << "fig_paper/fig_" << nucleus.c_str() << "_br.pdf";
	c->Print(os.str().c_str());

	if(nucleus=="11B"){
		TCanvas* c_br = new TCanvas("c_br","c_br",0,0,1200,600);
		c_br->Divide(3,2);
		int index=0;
		for(int i=0;i<part_num;i++){
			if(i==5) continue; //he
			double x,y;
			g[i]->GetPoint(bin,x,y);
			index++;
			c_br->cd(index);
			gPad->SetRightMargin(0.02);
			gPad->SetTopMargin(0.02);
			gPad->SetLeftMargin(0.13);
			gPad->SetBottomMargin(0.13);
			TH1F* waku_br = gPad->DrawFrame(0,0,25,0.5);
			waku_br->GetXaxis()->SetTitle("Excitation energy (MeV)");
			waku_br->GetXaxis()->CenterTitle();
			waku_br->GetYaxis()->SetTitle("Relative branching ratio");
			waku_br->GetYaxis()->CenterTitle();
			waku_br->GetYaxis()->SetNdivisions(505);
			g_ex[i]->Draw("L");
			TLatex* l = new TLatex(10,0.38,part_name[i]);
			l->SetTextFont(12);
			l->SetTextSize(0.08);
			l->Draw("same");
			os.str("");
			os << "Absolute Br:" << scientific << setprecision(1) << y;
			TLatex* l_br = new TLatex(10,0.30,os.str().c_str());
			//l_br->SetTextFont(12);
			l_br->SetTextSize(0.08);
			l_br->Draw("same");
			//cout << y << endl;
		}
		os.str("");
		os << "fig_paper/fig_" << nucleus.c_str() << "_br_Ex.pdf";
		c_br->Print(os.str().c_str());
	}

	return 0;
}
