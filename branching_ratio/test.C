int test(){

	TGraph* g;
	if(g==NULL) cout << "NULL" << endl;
	cout << g << endl;
	TFile* rootf = new TFile("output/Br_11B.root","READ");
	g = (TGraph*) rootf->Get("g_11B_br_0");
	cout << g << endl;

	rootf->Close();
	delete rootf;
	cout << g << endl;
	g->Draw();
	
	return 0;
}
