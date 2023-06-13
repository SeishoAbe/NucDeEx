#include "Particle.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <TLorentzVector.h>

///////////////
Particle::Particle()
///////////////
{
	_PDG=0;
	_mass=0;
	_momentum.SetXYZ(0,0,0);
	verbose=0;
}

///////////////
Particle::Particle(const int PDG, const double mass, const TVector3& mom, int v)
///////////////
{
	_PDG=PDG;
	_mass=mass;
	_momentum.SetXYZ(mom.X(), mom.Y(), mom.Z());
	verbose=v;
}


///////////////
double Particle::kE()
///////////////
{
	return sqrt( pow(_momentum.Mag(),2) + pow(_mass,2) ) - _mass;
}

///////////////
double Particle::totalE()
///////////////
{
	return sqrt( pow(_momentum.Mag(),2) + pow(_mass,2) );
}

///////////////
void Particle::Boost(const double totalE_parent,const TVector3& mom_parent)
///////////////
{ 
	// CM frame: parent nucleus 
	// totalE_parent: total energy of parent nucleus (kinetic E + mass E, but w/o excitation E)
	// mom_parent: momentum of parent nucleus:

	if(mom_parent.Mag()==0) return; // No boost (CM is at rest). Nothing to do

	// beta = v/c = momentum / total energy
	TVector3 beta;
	beta =1./totalE_parent*mom_parent;

	TLorentzVector lv(_momentum,kE()+_mass); // (TVector3, totalE)
	if(_mass<0 ||
			(_mass>0 && (lv.M()-_mass)/_mass > check_criteria)){
		cerr << "Unexpected TLorentzVector behaviour" << endl;
		cerr << "Mass from TLorentzVector = " << lv.M() <<endl;
		cerr << "Inputed Mass = " << _mass << endl;
		abort();
	}
	if(verbose>0) cout << "Bef boost: "; lv.Print();
	lv.Boost(beta);
	if(verbose>0) cout << "Aft boost: "; lv.Print();
	_momentum = lv.Vect();
	if(verbose>0) cout << "Aft boost:"; _momentum.Print();
	//lv.SetPxPyPzE(mom_parent.X(),mom_parent.Y(),mom_parent.Z(),total_E_parent);
}
