#ifndef __PARTICLE__HH__
#define __PARTICLE__HH__

#include "consts.hh"
#include <TVector3.h>

using namespace std;

class Particle{
	public:
	Particle();
	Particle(const int PDG,const double mass, const TVector3& mom,
					 const string name, const bool flag=1, const double Ex=0,
					 const int v=1);
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
	bool _flag;
	//	1 -> Needs to be tracked
	//  0 -> Please ignore it. It's intermediate state
	double _Ex;
	// This is only for intermediate state.
	// This shoud be zero if _flag==1.

	private:
	const double check_criteria=1e-3;
	int verbose;
};
#endif
