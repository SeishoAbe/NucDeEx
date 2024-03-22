#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "NucDeExUtils.hh"
#include "NucDeExEventInfo.hh"

/////////////////////////////////////////////
NucDeExEventInfo::NucDeExEventInfo(): eventID(0), fStatus(0), fShell(0), Zt(0), Nt(0), Z(0), N(0),
                                      Ex(0), ParticleVector(0)
/////////////////////////////////////////////
{
  Pinit.SetXYZ(0,0,0);
}

/////////////////////////////////////////////
NucDeExEventInfo::~NucDeExEventInfo()
/////////////////////////////////////////////
{
  if(ParticleVector!=0){
    ParticleVector->clear();
    delete ParticleVector;
  }
}
/////////////////////////////////////////////
void NucDeExEventInfo::InitParameters()
/////////////////////////////////////////////
{
  eventID = fStatus = fShell =0;
  Zt = Nt = Z = N = 0;
  Ex = 0.;
  if(ParticleVector!=0){
    ParticleVector->clear();
    delete ParticleVector;
  }
  ParticleVector = new std::vector<NucDeExParticle>;
}
