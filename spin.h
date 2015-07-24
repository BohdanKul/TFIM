#ifndef SPIN_H
#define SPIN_H

#include "randombase.h"
#include <iostream>
#include <vector>

using namespace std;
/******************************************************************************
 * Base spins class. 
 * Defines spins vector structure and basic operations on it.
 *****************************************************************************/
class Spins
{
    public:
        Spins(int _Nspins, long _seed); // initialize spins to a random distribution of +1/-1 
        Spins(int _Nspins, long _seed, vector<int>& _clamped); // set clamped spins to preset values 
        int  getSpin(int index);           // get the spin state
        int  getSize() 
            {return Nspins;};              // get the total number of spins 
        bool isClamped(int index)          // check whether a spin is clamped 
            {return bspins[index];};
        void setSpin(int index, int val);  // set the spin state
        void flip(int index);              // flip spin
        void print();                      // print out the spins state
    
        //Spins operator=(const Spins& clone)
        //    {Nsites = clone.Nsites; spins = clone.spins; return *this;};

    protected:
        vector<int>   spins;  // current spins state
        vector<bool> bspins;  // clamped or not-clamped state  
        int Nspins;           // number of spins
};


#endif
