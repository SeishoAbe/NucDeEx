#include <iomanip>

using namespace std;

int cout_energy(){
	const int bin=100;
	const double min=0, max=100;
	TH1D* h = new TH1D("h","",bin,min,max);

	ofstream ofs("energy.1.2.p.txt");
	
	ofs << setw(5) << bin << " 9 2" << endl;
	// (energy bin) (#spin) (#parity)


	// first, negative parity
	for(int i=0;i<bin;i++){
		double e = h->GetBinCenter(i+1);
		//ofs << setw(5) << setfill('0') << e 
		ofs << setw(5) << e 
        << " 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" << endl;
				// all zero popolation
	}
	// then, positive parity
	for(int i=0;i<bin;i++){
		double e = h->GetBinCenter(i+1);
		//ofs << setw(5) << setfill('0') << e 
		ofs << setw(5) << e 
        << " 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" << endl;
				// all zero popolation
	}
	ofs.close();

	return 0;
}
