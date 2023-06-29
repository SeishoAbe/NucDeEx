int plot_15N(){
	// scanned (Fig3)
	// data
	double br_n_2b = 9.8;
	double br_n    = 24.7;
	double br_p_2b = 4.5;
	double br_p    = 17.7;
	double br_d_2b = 6.9;
	double br_d    = 9.7;
	double br_t_2b = 6.8;
	double br_t    = 8.7;
	double br_a_2b = 2.4;
	double br_a    = 4.7;
	// error
	double br_e_n_2b = 1.0;
	double br_e_n    = 1.6;
	double br_e_p_2b = 0.6;
	double br_e_p    = 1.1;
	double br_e_d_2b = 0.4; 
	double br_e_d    = 0.5;
	double br_e_t_2b = 0.4;
	double br_e_t    = 0.4;
	double br_e_a_2b = 0.5;
	double br_e_a    = 0.8;

	// cascade prediction
	double cas_br_n_2b = 5.1;
	double cas_br_n    = 23.9;
	double cas_br_p_2b = 2.8;
	double cas_br_p    = 19.3;
	double cas_br_d_2b = 5.7;
	double cas_br_d    = 11.5;
	double cas_br_t_2b = 1.4;
	double cas_br_t    = 2.4;
	double cas_br_a_2b = 3.0;
	double cas_br_a    = 7.2;

	TFile* rootf = new TFile("Yosoi_15N.root","RECREATE");
	TEnv env("env");
	env.SetValue("br_p_2b",br_p_2b);
	env.SetValue("br_p",br_p);
	env.SetValue("br_d_2b",br_d_2b);
	env.SetValue("br_d",br_d);
	env.SetValue("br_t_2b",br_t_2b);
	env.SetValue("br_t",br_t);
	env.SetValue("br_a_2b",br_a_2b);
	env.SetValue("br_a",br_a);
	//
	env.SetValue("br_e_p_2b",br_e_p_2b);
	env.SetValue("br_e_p",br_e_p);
	env.SetValue("br_e_d_2b",br_e_d_2b);
	env.SetValue("br_e_d",br_e_d);
	env.SetValue("br_e_t_2b",br_e_t_2b);
	env.SetValue("br_e_t",br_e_t);
	env.SetValue("br_e_a_2b",br_e_a_2b);
	env.SetValue("br_e_a",br_e_a);
	//
	env.SetValue("cas_br_n_2b",cas_br_n_2b);
	env.SetValue("cas_br_n",cas_br_n);
	env.SetValue("cas_br_p_2b",cas_br_p_2b);
	env.SetValue("cas_br_p",cas_br_p);
	env.SetValue("cas_br_d_2b",cas_br_d_2b);
	env.SetValue("cas_br_d",cas_br_d);
	env.SetValue("cas_br_t_2b",cas_br_t_2b);
	env.SetValue("cas_br_t",cas_br_t);
	env.SetValue("cas_br_a_2b",cas_br_a_2b);
	env.SetValue("cas_br_a",cas_br_a);
	env.Print();
	env.Write("env");


	return 0;
}
