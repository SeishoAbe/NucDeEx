int plot(){
	const double rbr_nda_n=71.4462;
	const double rbr_nda_da=28.5538;


	const double br_n_2b = 9.33421;
	const double br_n    = 21.9994;
	const double br_p_2b = 5.70191;
	const double br_p    = 7.09744;
	const double br_d_2b = 4.11364;
	const double br_d    = 5.67756;
	const double br_t_2b = 1.51171;
	const double br_t    = 2.17326;
	const double br_a_2b = 3.34729;
	const double br_a    = 3.42384;

	TFile* rootf = new TFile("KamLAND.root","RECREATE");
	TEnv env("env");
	env.SetValue("rbr_nda_n",rbr_nda_n);
	env.SetValue("rbr_nda_da",rbr_nda_da);
	//
	env.SetValue("br_n_2b",br_n_2b);
	env.SetValue("br_n",br_n);
	env.SetValue("br_p_2b",br_p_2b);
	env.SetValue("br_p",br_p);
	env.SetValue("br_d_2b",br_d_2b);
	env.SetValue("br_d",br_d);
	env.SetValue("br_t_2b",br_t_2b);
	env.SetValue("br_t",br_t);
	env.SetValue("br_a_2b",br_a_2b);
	env.SetValue("br_a",br_a);
	env.Print();
	env.Write("env");


	return 0;
	
}
