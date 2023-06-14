#include "include/consts.hh"
int plot_simulation(){
	TFile* rootf = new TFile("sim_out/output_11B.root","READ");
	TTree* tree = (TTree*) rootf->Get("tree");
	int eventID, size;
	double MissE, Ex, S;
	double PinitMag, PinitX,PinitY,PinitZ;
	//
	int PDG[bins];
	double mass[bins];
	double totalE[bins],kE[bins];
	double PMag[bins], PX[bins],PY[bins],PZ[bins];
	TBranch *b_PDG, *b_mass,*b_totalE, *b_kE;
	TBranch *b_PMag, *b_PX, *b_PY, *b_PZ;
	tree->SetBranchAddress("eventID",&eventID);
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
	
	for(int i=0;i<tree->GetEntries();i++){
		tree->GetEntry(i);
		cout << PMag[0] << endl;
	}
	
	tree->Draw("PMag[0]");

	return 0;
}
