#include "spins.h"
#include "lattice.h"

using namespace std;

class Hamiltonian
{
    public
        Hamiltonian(Spins* _spins, Bonds* _bonds, Field* _zfield, Field* _xfield);
        pair<int,int> getSpins(int index);   // get spins associated with a bond
        void flipBondSpins(int index);       // flip spins belonging to a bond
        void computeDiagProb();
        
        Spins& spins;
        Bonds& bonds;
        Field& zfield;
        Field& xfield; 

        float getEoffset(){     return tEoffset;};
        float getBondEoffset(){ return bEoffset;};
        float getEtotal(){      return tE; } 
    private
        float tEoffset; // energy off-set due to all diagonal operators
        float bEoffset; // energy off-set due to bond diagonal operators only
        float tE;       // total diagonal energy
}
