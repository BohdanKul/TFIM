#include "hamiltonian.h"
#include <cmath>

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
    sites.flipSiteSpin(sites.first);
    sites.flipSiteSpin(sites.second);
}

/*****************************************************************************
 * Initialize spins and interactions, 
 * compute energy offsets due to diagonal operators
*****************************************************************************/
float Hamiltonian::Hamiltonian():
spins(*_spins), bonds(*_bonds), zfield(*_zfield), xfield(*_xfield)
{
    // Compute the energy off-set due to diagonal operators
    tEoffset = 0;
    for (auto elem = bonds.begin(); elem!=bonds.end(); elem++)
        tEoffset += abs(elem->getStrength());
    bEoffset = tEoffset; 

    for (auto elem = zfield.begin(); zfield!=bonds.end(); elem++)
        tEoffset += *elem;

    tE = tEoffset + bEoffset;
}
