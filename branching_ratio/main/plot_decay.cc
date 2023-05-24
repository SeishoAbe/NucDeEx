#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 


#include "NucleusTable.hh"
#include "ReadTALYS.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>

using namespace std;

int main(int argc, char* argv[]){
	if(argc!=2){
		cerr << argv[0] << " [TALYS output]" << endl;
		return 0;
	}
  // --- FIXME  --- //
  std::string parent="11C";
  const int maxlevelstar=8;
  const int maxlevelsbin_g=8;
  const int maxlevelsbin_n=10;
  const int maxlevelsbin_p=10;
  const int maxlevelsbin_a=10;
  const int maxlevelsbin_d=5;
  const int maxlevelsbin_t=5;
  const int maxlevelsbin_h=5;

  const int Ex_max=60000; // keV
  const int Ex_bin_width=100; // keV
  const double max_Ex_plot=50;
  // ---------------//
  
  std::ostringstream os;
  NucleusTable* nucleus_table = new NucleusTable();
  if(!nucleus_table->ReadTables()){
		cerr << "something wrong" << endl;
		return 1;
	}
	
	//ReadTALYS* read_talys = new ReadTALYS();//ReadTALYS(argv[1], nucleus_table);
	ReadTALYS* read_talys = new ReadTALYS(argv[1], nucleus_table);
	read_talys->Read();


	delete read_talys;
	delete nucleus_table;
  return 0;
}
