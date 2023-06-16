int plot(){

	// --- Measured data (Table 1) --- //
	const double br_p_2b = 8.5;
	const double br_p    = 11.7;
	const double br_d_2b = 5.3;
	const double br_d    = 6.5;
	const double br_t_2b = 16.4;
	const double br_t    = 18.6;
	const double br_a_2b = 7.7;
	const double br_a    = 16.8;
	// error
	const double br_e_p_2b = 0.5;
	const double br_e_p    = 0.6;
	const double br_e_d_2b = 0.3;
	const double br_e_d    = 0.4;
	const double br_e_t_2b = 0.6;
	const double br_e_t    = 0.6;
	const double br_e_a_2b = 0.4;
	const double br_e_a    = 0.6;

	// --- CASCADE prediction (scanned) --- //
	const double cas_br_n_2b = 12.0;
	const double cas_br_n    = 21.4;
	const double cas_br_p_2b = 5.5;
	const double cas_br_p    = 9.9;
	const double cas_br_d_2b = 8.2;
	const double cas_br_d    = 10.5;
	const double cas_br_t_2b = 4.1;
	const double cas_br_t    = 4.5;
	const double cas_br_a_2b = 5.0;
	const double cas_br_a    = 8.1;

	TFile* rootf = new TFile("Yosoi.root","RECREATE");
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
