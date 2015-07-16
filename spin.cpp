#include "spin.h"


/*****************************************************************************
 * Initializes spins in a random state without underlying geometry
 *****************************************************************************/
Spins::Spins(int _Nspins, long _seed)
/*
    _Nspins - sets the size of spins vectors
    _seed   - initializes the random number generator necessary 
              to generate a random distribution of spins
*/
{
    RandomBase  rand(_seed);     // random numbers generator
   
    // Initialize spins randomly with values 1 and -1
    Nspins = _Nspins;
    spins.resize(Nspins,0);
    for (auto spin=spins.begin(); spin!=spins.end(); spin++) {                          
        *spin = pow(-1,rand.uRandInt()%2);                         
    }

    cout << endl << "=== Spins initialization ===" << endl;
    cout << "   spins #: " << Nspins << endl;
    cout << "   state  : ";
    for (auto spin=spins.begin(); spin!=spins.end(); spin++) {                          
        cout << *spin << " "; 
    }
    cout << endl << endl;

}


/*****************************************************************************
 * Returns the spin state 
 *****************************************************************************/
int Spins::getSpin(int index)
/*
   index - index of the element to be flipped  
*/
{
    return spins[index];
}

/*****************************************************************************
 * Sets the spin state to a particular value
 *****************************************************************************/
void Spins::setSpin(int index, int val)
/*
   index - index of the element to be flipped  
   val   - the new value of the spin state
*/
{
    spins[index] = val;
}

/*****************************************************************************
 * Flips a spin
 *****************************************************************************/
void Spins::flip(int index)
/*
   index - index of the element to be flipped  
*/
{
    setSpin(index, -1*getSpin(index));
}

/*****************************************************************************
 * Print out the spins state 
 *****************************************************************************/
void Spins::print()
{
    cout << "   Spins state: " ;
    for (int i=0; i!=Nspins; i++)
        cout << getSpin(i) << " ";
    cout << endl;
}

//*****************************************************************************
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//*****************************************************************************


