#ifndef __NUCDEEXRANDOM__HH__
#define __NUCDEEXRANDOM__HH__

#include <TRandom3.h>

namespace NucDeEx{
  namespace Random{
    extern TRandom3 rndm;
    // needs "extern". These are defined in *.cc

    double random();
    void SetSeed(int s);
  }
}
#endif
