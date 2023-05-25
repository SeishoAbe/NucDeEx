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
#endif
