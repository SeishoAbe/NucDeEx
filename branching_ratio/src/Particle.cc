#include "Particle.hh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cstdlib>

///////////////
Particle::Particle()
///////////////
{
	_PDG=0;
	_mass=0;
	_momentum = new TVector3(0,0,0);
}

///////////////
Particle::Particle(int PDG,double mass, TVector3* momentum)
///////////////
{
	_PDG=PDG;
	_mass=mass;
	_momentum = new TVector3(0,0,0);
	_momentum->SetXYZ(momentum->X(), momentum->Y(), momentum->Z());
}

///////////////
Particle::~Particle()
///////////////
{
	delete _momentum;
}


///////////////
double Particle::KineticE()
///////////////
{
	return sqrt( pow(_momentum->Mag(),2) + pow(_mass,2) ) - _mass;
}

///////////////
void Particle::Boost()
///////////////
{
	;
}
