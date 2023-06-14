#ifndef __CONSTS__HH__
#define __CONSTS__HH__

#include <string>

using namespace std;

static const int bins=100;
static const int parity=2;

static const int num_particle=7;
static const string particle_name[num_particle]
	= {"gamma","neutron","proton",
		 "deuteron", "triton","helium-3","alpha"};

static const int PDG_particle[num_particle]
	= {22, 2112, 2212,
		 1000010020, 1000010030, 1000020030,1000020040};
static const int color_root[num_particle]
	= {880+1, // violet+1
		 600, // bule
		 632, // red
		 920+1, // gray +1
		 432+1, // cyan+1
		 416+1,// green+1
		 616 // magenta
		};

#endif
