#include "hamiltonian.h"
#include <cmath>


/*****************************************************************************
 * Initialize spins and interactions, 
 * compute energy offsets due to diagonal operators
*****************************************************************************/
Hamiltonian::Hamiltonian(Spins* _spins, Bonds* _bonds, vector<float>* _xfield):
spins(*_spins), bonds(*_bonds), xfield(*_xfield)
{
    cout << "---Hamiltonian initialization---" << endl;
    cout << "   Transverse field:" << endl << "   ";
    for (auto xf=xfield.begin(); xf!=xfield.end(); xf++){
        cout << *xf << " ";
    }
    cout << endl;

    // Compute the energy off-set due to diagonal operators
    tEoffset = 0;
    for (int index = 0; index!=bonds.getBondsN(); index++)
        tEoffset += abs(bonds.getStrength(index));
    bEoffset = tEoffset; 

    for (auto elem = xfield.begin(); elem!=xfield.end(); elem++)
        tEoffset += *elem;

    tE = tEoffset + bEoffset;
    cout << "   Bond offset: "           << getBondEoffset() << endl 
         << "   Total offset: "          << getEoffset()     << endl
         << "   Total diagonal energy: " << getEtotal()      << endl;
    cout << endl;
}

/*****************************************************************************
 * Return the state of two spins associated with a bond
 *****************************************************************************/
pair<int,int> Hamiltonian::getSpins(int index)
/*
    index - index of the bond 
*/
{
    pair<int,int> sites = bonds.getSites(index);
    return pair<int,int>  (spins.getSpin(sites.first),spins.getSpin(sites.second));  
}


/*****************************************************************************
 * Flip both spins associated with a bond 
 *****************************************************************************/
void Hamiltonian::flipBondSpins(int index)
/*
    index - index of the bond 
*/
{
    pair<int,int> sites = bonds.getSites(index);
    spins.flip(sites.first);
    spins.flip(sites.second);
}


