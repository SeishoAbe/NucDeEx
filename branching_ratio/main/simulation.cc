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
	double Ex=23; // MeV
	int Z=5;
	int N=6;

	Deexcitation* deex = new Deexcitation();
	//deex->SetSeed(0); // 0 -> time
	deex->SetSeed(1);
	deex->SetVerbose(1);

	for(int i=0;i<10;i++){
		deex->DoDeex(Z,N,Ex);
	}


	return 0;
}