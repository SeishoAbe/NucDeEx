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

#include "NucDeExEventInfo.hh"

/////////////////////////////////////////////
NucDeExEventInfo::NucDeExEventInfo(): EventID(0), fStatus(0), fShell(0),
                                      Zt(0), Nt(0), Z(0), N(0), Ex(0)
/////////////////////////////////////////////
{
  Pinit.SetXYZ(0,0,0);
  ParticleVector.clear();
}

/////////////////////////////////////////////
NucDeExEventInfo::~NucDeExEventInfo()
/////////////////////////////////////////////
{
  ;
}
/////////////////////////////////////////////
void NucDeExEventInfo::InitParameters()
/////////////////////////////////////////////
{
  EventID = fStatus = fShell =0;
  Zt = Nt = Z = N = 0;
  Ex = 0.;
  ParticleVector.clear();
  std::vector<NucDeExParticle>().swap(ParticleVector);
}
