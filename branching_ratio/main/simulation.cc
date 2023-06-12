#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 

#include "Deexcitation.hh"

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
	// paremater 
	double Ex=22; // MeV
	int Z=5;
	int N=6;

	Deexcitation* deex = new Deexcitation();
	deex->SetSeed(0); // 0 -> time
	deex->SetVerbose(1);
	deex->DoDeex(Z,N,Ex);


	return 0;
}
