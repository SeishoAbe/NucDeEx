#ifndef __CONSTS__HH__
#define __CONSTS__HH__

#include <string>

using namespace std;

static const int bins=100;
static const int parity=2;
static const char* nuc_name[]={"","H","He","Li","Be","B","C","N","O"};
	// [Z]

static const int num_particle=7;
static const string particle_name[num_particle]
	= {"gamma","neutron","proton",
		 "deuteron", "triton","helium-3","alpha"};
static const int color_root[num_particle]
	= {880+1, // violet+1
		 600, // bule
		 632, // red
		 920+1, // gray +1
		 432+1, // cyan+1
		 416+1,// green+1
		 616 // magenta
		};

static const string decay_name[num_particle] // for G4
	= {"IT","Neutron","Proton",
		 "Deuteron","Triton","He3","Alpha"};
static const float check_criteria=0.05;


//static const float Ex_max=100*1e3; // keV
//static const float Ex_bin_width=0.2*1e3; // keV
#endif
