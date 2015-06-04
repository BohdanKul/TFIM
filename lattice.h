#ifndef LATTICE_H
#define LATTICE_H

#include <string.h>
#include <vector>
#include "randombase.cpp"
using namespace std;


/******************************************************************************
 * Base spins class. 
 * Defines spins vector structure and basic operations on it.
 *****************************************************************************/
class Spins
{
    public:
        Spins(int _Nsites);
        Spins(int _Nsites, long _seed);
        //Spins operator=(const Spins& clone)
        //    {Nsites = clone.Nsites; spins = clone.spins; return *this;};
        int  getSpin(int index);           // get the spin state
        void setSpin(int index, int val);  // set the spin value
        void flipSiteSpin(int index);      // flip a particular spin
        int  getSize() {return Nsites;};
        void printSpins();
    //protected:
        int Nsites;           // number of sites
        vector<int> spins;    // list of spin variables
};



/******************************************************************************
 *A container class holding information about 2 sites of a bond
 *****************************************************************************/
class tBond
{
    public:
        int A;
        int B;
        tBond(): A(-1), B(-1) {};
        tBond(int x, int y): A(x), B(y) {};
        tBond set(int x, int y){ A = x; B = y; return *this;};
};

/******************************************************************************
 * Base lattice class. 
 * Defines bonds operating on top of the spins structure it inherites. 
 *****************************************************************************/
class Lattice: public Spins
{
    public:
        Lattice(int _x, int _y, int _unitx, int _unity, long _seed);
        pair<int,int> getSites(int index);   // get sites associated with a bond
        pair<int,int> getSpins(int index);   // get spins associated with a bond
        void flipBondSpins(int index);       // flip spins belonging to a bond
        void setBond(int siteA, int siteB);  // create a bond between two sites
        void printBonds();
        int  getBondsN(){ return Nbonds;};
        string& getName(){return name;};

        int unitx;  // unit cell width
        int unity;  // unit cell height
        int x;      // lattice width (in terms of unit cells)
        int y;      // lattice height
    protected:
        string name;
        vector<tBond> bonds;  // list of bonds between those spins  
        int Nbonds; // number of bonds

};


/******************************************************************************
 * Rectangular lattice class
 *****************************************************************************/
class Rectangle: public Lattice
{
    public:
        Rectangle(int _x, int _y, bool _OBC, long _seed);
};

/******************************************************************************
 * Chimera lattice class
 *****************************************************************************/
class Chimera: public Lattice
{
    public:
        Chimera(int _x, int _y, long _seed,  int _unity);
};

#endif
