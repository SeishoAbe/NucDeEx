#ifndef __NUCDEEXRANDOM__HH__
#define __NUCDEEXRANDOM__HH__

#include <TRandom3.h>

namespace NucDeEx{
  namespace Random{
    TRandom3* rndm = NULL;
    double random();
  }
}

#endif
