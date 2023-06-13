#ifndef __PARTICLE__HH__
#define __PARTICLE__HH__

#include "consts.hh"
#include <TVector3.h>

using namespace std;

class Particle{
	public:
	Particle();
	Particle(int PDG,double mass, TVector3* momentum);
	~Particle();
	
	//double momentum(){_momentum->Mag();};
	double KineticE();

	void Boost();

	//private:
	int _PDG;
	double _mass;
	TVector3* _momentum;
};
#endif
