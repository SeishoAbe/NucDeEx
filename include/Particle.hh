#ifndef __PARTICLE__HH__
#define __PARTICLE__HH__

#include "consts.hh"
#include <TVector3.h>

using namespace std;

class Particle{
	public:
	Particle();
	Particle(const int PDG,const double mass, const TVector3& mom, const string name);
		// TVector3 is called by referecnce
	~Particle(){;};
	
	double kE(); //kinetic energy
	double totalE(); // total energy

	void Boost(const double totalE_parent, const TVector3& mom_parent);

	//private:
	int _PDG;
	double _mass;
	TVector3 _momentum;
	string _name;

	private:
	const double check_criteria=1e-3;
	int verbose;
};
#endif
