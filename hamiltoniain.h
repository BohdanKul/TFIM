#include "spins.h"
#include "lattice.h"

using namespace std;

class Hamiltonian
{
    public
        Hamiltonian(Spins* _spins, Bonds* _bonds, Field* _zfield, Field* _xfield):
                spins(*_spins), bonds(*_bonds), zfield(*_zfield), xfield(*_xfield){};

        pair<int,int> getSpins(int index);   // get spins associated with a bond
        void flipBondSpins(int index);       // flip spins belonging to a bond
        

        Spins& spins;
        Bonds& bonds;
        Field& zfield;
        Field& xfield;  
}
