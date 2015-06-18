#include "randombase.h"
using namespace std;

RandomBase::RandomBase(long seed):
eng(1),uReal(0,1),uInt(0,655352310000),uRandInt(eng,uInt),uRand(eng,uReal)
{    
    //Initialize ranndomness with the seed
    uRand.engine().seed(seed);
    uRand.distribution().reset();
    uRandInt.engine().seed(seed);
    uRandInt.distribution().reset();
}

