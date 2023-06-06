#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc,char* argv[]){
	if(argc!=2){
		cerr << "Input: " << argv[0] << " [Target nuclues]" << endl;
		return 1;
	}

	string keyword="):";
	ostringstream os;
	os.str("");
	os << "output/output_" << argv[1];
	ifstream ifs(os.str().c_str());
	if(!ifs.is_open()){
		cerr << "Cannot read " << os.str().c_str() << endl;
		return 1;
	}
	
	os.str("");
	os << "separation_energy_" << argv[1] << ".txt";
	ofstream ofs(os.str().c_str());

	char buf[500];
	while(ifs.getline(buf,sizeof(buf))){
		string st = string(buf);
		int find = st.find(keyword.c_str()); // Q(X,X): 
		if(find!=string::npos){
			st = st.substr(find+keyword.length());
			float S;
			istringstream(st) >> S;
			S *= -1; // (usually) negative -> positive (just a notation..)
			if(S==0) S=0;
			//cout << S << endl;
			ofs << S << endl;
		}
	}
	ofs.close();

	
	return 0;
}
