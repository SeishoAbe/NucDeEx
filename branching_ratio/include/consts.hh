#ifndef __CONSTS__HH__
#define __CONSTS__HH__

#include <string>

using namespace std;

static const int bins=100;
static const int parity=2;
static const int num_particle=7;
static const string particle_name[num_particle]
	= {"gamma","neutron","proton","alpha",
		 "deuteron", "triton","helium-3"};
static const int color_root[num_particle]
	= {880+1, // violet+1
		 600, // bule
		 632, // red
		 616, // magenta
		 920+1, // gray +1
		 432+1, // cyan+1
		 416+1 // green+1
		};

static const float check_criteria=0.02;
#endif
