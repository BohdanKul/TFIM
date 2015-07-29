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
        int  getSpin(int index);           // get the spin state
        int  getSize() {return Nspins;};   // get the total number of spins 
        void setSpin(int index, int val);  // set the spin state
        void flip(int index);              // flip spin
        void print();                      // print out the spins state
        vector<int>* getState(){ return &spins;}; 
        //Spins operator=(const Spins& clone)
        //    {Nsites = clone.Nsites; spins = clone.spins; return *this;};

    protected:
        vector<int> spins;    // current spins state
        int Nspins;           // number of spins
};


#endif
