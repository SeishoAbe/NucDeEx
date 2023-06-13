#ifndef __PARTICLE__HH__
#define __PARTICLE__HH__

#include "consts.hh"
#include <TVector3.h>

using namespace std;

class Particle{
	public:
	Particle();
	Particle(const int PDG,const double mass, const TVector3& mom);
		// TVector3 is called by referecnce
	~Particle(){;};
	
	//double momentum(){_momentum->Mag();};
	double KineticE();

	void Boost();

	//private:
	int _PDG;
	double _mass;
	TVector3 _momentum;
};
#endif
