#include "../include/consts.hh"

using namespace std;
int cout_br_for_hu(){
	// ---- FIXME ---- //
	const int ldmodel=3;
	const bool parity_optmodall=1;
	// --------------- //
	ostringstream os;
	os.str("");
	os << "sim_out/11B_ldmodel" << ldmodel;
	if(parity_optmodall) os << "_parity_optmodall";
	os << "_Exmin15.9_raw.txt";
	cout << os.str().c_str() << endl;
	ifstream ifs(os.str().c_str());
	if(!ifs.is_open()) return 1;


	char buf[256];
	
	const int nmode=5+1;
	const string mode_hu[nmode]
		= {"n10B","p10Be","d9Be","t8Be","a7Li","nn9B"};
	double br[nmode]={0};

	double br_np9Be=0, br_nd8Be=0, br_na6Li=0, br_da5He=0, br_ta4He=0;
	double br_nnp8Be=0, br_npp8Li=0, br_nnd7Be=0;
	double br_npd7Li=0, br_npt6Li=0;
	double br_nna5Li=0, br_npa5He=0;
	double br_nnpp7Li=0, br_nnpd6Li=0;
	double br_others=0;
	while(ifs.getline(buf,sizeof(buf))){
		string st = buf;
		string mode;
		double br_tmp;
		istringstream(buf) >> mode >> br_tmp;
		//
		bool flag=0;
		for(int i=0;i<nmode;i++){
			if(mode==mode_hu[i]){
				br[i]=br_tmp;
				flag=1;
				break;
			}
		}
		if(flag) continue;
		if(mode=="np9Be" || mode=="pn9Be") br_np9Be+=br_tmp;
		else if(mode=="nd8Be" || mode=="dn8Be") br_nd8Be+=br_tmp;
		else if(mode=="na6Li" || mode=="an6Li") br_na6Li+=br_tmp;
		else if(mode=="da5He" || mode=="ad5He") br_da5He+=br_tmp;
		else if(mode=="ta4He" || mode=="at4He") br_ta4He+=br_tmp;

		else if(mode=="nnp8Be" || mode=="npn8Be" || mode=="pnn8Be") br_nnp8Be+=br_tmp;
		else if(mode=="npp8Li" || mode=="pnp8Li" || mode=="ppn8Li") br_npp8Li+=br_tmp;
		else if(mode=="nnd7Be" || mode=="ndn7Be" || mode=="dnn7Be") br_nnd7Be+=br_tmp;

		else if(mode=="npd7Li" || mode=="ndp7Li" || mode=="pnd7Li"
				|| mode=="pdn7Li" || mode=="dnp7Li" || mode=="dpn7Li") br_npd7Li+=br_tmp;
		else if(mode=="npt6Li" || mode=="ntp6Li" || mode=="pnt6Li"
				|| mode=="ptn6Li" || mode=="tnp6Li" || mode=="tpn6Li") br_npt6Li+=br_tmp;

		else if(mode=="nna5Li" || mode=="nan5Li" || mode=="ann5Li") br_nna5Li+=br_tmp;
		else if(mode=="npa5He" || mode=="nap5He" || mode=="pna5He"
				|| mode=="pan5He" || mode=="anp5He" || mode=="apn5He") br_npa5He+=br_tmp;

		else if(mode=="nnpp7Li" || mode=="npnp7Li" || mode=="nppn7Li"
				|| mode=="ppnn7Li" || mode=="pnpn7Li" || mode=="pnnp7Li") br_nnpp7Li+=br_tmp;
		else if(mode=="nnpd6Li" || mode=="nndp6Li" || mode=="ndnp6Li"  
				|| mode=="ndpn6Li" || mode=="npdn6Li" || mode=="npnd6Li"
				|| mode=="pnnd6Li" || mode=="pndn6Li" || mode=="pdnn6Li"
				|| mode=="dnnp6Li" || mode=="dpnn6Li" || mode=="dnpn6Li") br_nnpd6Li+=br_tmp;
		else br_others+=br_tmp;
	}


	for(int i=0;i<nmode;i++){
		cout << mode_hu[i] << " " << br[i] << endl;
	}
	cout << "np9Be  " << br_np9Be << endl;
	cout << "nd8Be  " << br_nd8Be << endl;
	cout << "na6Li  " << br_na6Li << endl;
	cout << "da5He  " << br_da5He << endl;
	cout << "ta4He  " << br_ta4He << endl;
	//
	cout << "nnp8Be  " << br_nnp8Be << endl;
	cout << "npp8Li  " << br_npp8Li << endl;
	cout << "nnd7Be  " << br_nnd7Be << endl;
	//
	cout << "npd7Li  " << br_npd7Li << endl;
	cout << "npt6Li  " << br_npt6Li << endl;
	//
	cout << "nna5Li  " << br_nna5Li << endl;
	cout << "npa5He  " << br_npa5He << endl;
	//
	cout << "nnpp7Li  " << br_nnpp7Li << endl;
	cout << "nnpd6Li  " << br_nnpd6Li << endl;
	cout << "others   " << br_others << endl;

	//
	for(int i=0;i<nmode;i++){
		cout << br[i] << endl;
	}
	cout << br_np9Be << endl;
	cout << br_nd8Be << endl;
	cout << br_na6Li << endl;
	cout << br_da5He << endl;
	cout << br_ta4He << endl;
	//
	cout << br_nnp8Be << endl;
	cout << br_npp8Li << endl;
	cout << br_nnd7Be << endl;
	//
	cout << br_npd7Li << endl;
	cout << br_npt6Li << endl;
	//
	cout << br_nna5Li << endl;
	cout << br_npa5He << endl;
	//
	cout << br_nnpp7Li << endl;
	cout << br_nnpd6Li << endl;
	cout << br_others << endl;

	return 0;
}
