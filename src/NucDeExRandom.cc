#include "NucDeExRandom.hh"

namespace NucDeEx{
  namespace Random{
    TRandom3 rndm(1);
    double random(){ return rndm.Rndm(); };
    void SetSeed(int s){ rndm.SetSeed(s); };
  }
}
