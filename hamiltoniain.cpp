#include "hamiltonian.h"

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
