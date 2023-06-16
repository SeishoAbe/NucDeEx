int plot(){
	const double rbr_nda_n=80.5;
	const double rbr_nda_da=19.5;

	const double br_n_2b = 8.8;
	const double br_n    = 25.5;
	const double br_p_2b = 5.0;
	const double br_p    = 14.4;
	const double br_d_2b = 5.9;
	const double br_d    = 10.8;
	const double br_t_2b = 4.0;
	const double br_t    = 4.4;
	const double br_a_2b = 4.6;
	const double br_a    = 16.2;

	TFile* rootf = new TFile("Hu.root","RECREATE");
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
