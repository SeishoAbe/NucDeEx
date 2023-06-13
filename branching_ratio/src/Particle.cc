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
	_momentum.SetXYZ(0,0,0);
}

///////////////
Particle::Particle(const int PDG, const double mass, const TVector3& mom)
///////////////
{
	_PDG=PDG;
	_mass=mass;
	_momentum.SetXYZ(mom.X(), mom.Y(), mom.Z());
}


///////////////
double Particle::KineticE()
///////////////
{
	return sqrt( pow(_momentum.Mag(),2) + pow(_mass,2) ) - _mass;
}

///////////////
void Particle::Boost()
///////////////
{
	;
}
