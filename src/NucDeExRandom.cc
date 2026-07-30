#include <iostream>
#include "NucDeExRandom.hh"

namespace NucDeEx{
  namespace Random{
    TRandom3* rndm = nullptr;
    double random(){ 
      if(!rndm) rndm = new TRandom3(1);
      return rndm->Rndm();
    }
    void SetSeed(int s)
    { 
      if(!rndm) rndm = new TRandom3(1);
      rndm->SetSeed(s);
      std::cout << "NucDeEx::Random::SetSeed(" << s << ")" << std::endl;
    }
  }
}
