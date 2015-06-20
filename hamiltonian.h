#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "spin.h"
#include "bond.h"

using namespace std;

class Hamiltonian
{
    public:
        Hamiltonian(Spins* _spins, Bonds* _bonds, vector<float>* _xfield);
        pair<int,int> getSpins(int index);   // get spins associated with a bond
        void flipBondSpins(int index);       // flip spins belonging to a bond
        void computeDiagProb();
        
        Spins& spins;
        Bonds& bonds;
        vector<float>& xfield; 

        float getEoffset(){     return tEoffset;};
        float getBondEoffset(){ return bEoffset;};
        float getEtotal(){      return tE; } 
    private:
        float tEoffset; // energy off-set due to all diagonal operators
        float bEoffset; // energy off-set due to bond diagonal operators only
        float tE;       // total diagonal energy
};

#endif
