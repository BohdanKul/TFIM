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
        Spins(int _Nspins);                // initialize spins to 1
        Spins(int _Nspins, long _seed);    // initialize spins to a random
                                           // distribution of +1/-1 
        int  getSpin(int index);           // get the spin state
        void setSpin(int index, int val);  // set the spin state
        int  getSize() {return Nspins;};   // get the total number of spins 
        void flipSiteSpin(int index);      // flip the spin
        void printSpins();                 // print out the spins state
    
        //Spins operator=(const Spins& clone)
        //    {Nsites = clone.Nsites; spins = clone.spins; return *this;};
    protected:
        int Nspins;           // number of sites
        vector<int> spins;    // current spins state
};



/******************************************************************************
 * Class holding coordinates of 2 sites belonging to a bond, 
 * as well as its strength value.
 *****************************************************************************/
class tBond
{
    public:

        tBond(): A(-1), B(-1), strength(-1); {};                  // initiate the bond to default values
        tBond(int x, int y): A(x), B(y), strength(-1) {};         // initiate bond's spins to user-defined indices  
        tBond(int x, int y, float s): A(x), B(y), strength(s) {}; // initiate bond's spins and bond's strength to 
                                                                  // user-defined values

        tBond set(int x, int y)            // reset bond's spins to new values
            { A = x; B = y; return *this;};
        tBond set(int x, int y, float s)   // reset bond's spins and value it to new values
            { A = x; B = y; strength = s; return *this;};
        pair<int,int> get()                // get the indices of bond's spins 
            {return pair<int,int>(A,B);};
        float  getStrength()               // get bond's strength
            {return strength;}; 
        void print()                       // print the bond
            {cout<<"("<<A<<","<<B<<")="<<strength<<"\n";};  
    protected:
        int A;          // indices of two sites belonging to the bond 
        int B; 
        float strength; // bond strength as dictated by underlying Hamiltonian
};

/******************************************************************************
 * Virtual lattice class that adds a bonds structure operating on 
 * top of the spins structure it inherites. This new structure 
 * defines physical bonds existing between spins. This class only
 * implements the necessary functionality to deal with this new
 * structure without actually setting it up. In order to define 
 * a physical lattice with real bonds, one needs to inherit this 
 * class and define them in the constructors of the new class. 
 * (see rectangle and chimera classes) 
 *****************************************************************************/
class Lattice: public Spins
{
    public:
        Lattice(int _x, int _y, int _unitx,  // initialize lattice geometry and
                int _unity, long _seed);     // set its spins state to a random state

        pair<int,int> getSites(int index);   // get sites associated with a bond
        pair<int,int> getSpins(int index);   // get spins associated with a bond
        void flipBondSpins(int index);       // flip spins belonging to a bond
        void setBond(int siteA, int siteB);  // create a bond between two sites
        
        void printBonds();                   // print out all bonds 
        int  getBondsN(){ return Nbonds;};   // get the number of bonds
        string& getName(){return name;};     // get the lattice name

        int getWidth(){return x;};           // get lattice parameters
        int getHeight(){return y;};          //
        int getUnitWidth(){return unitx;};   //
        int getUnitHeight(){return unity;};  //
    protected:
        string name;          // lattice-type name
        
        vector<tBond> bonds;  // vector of bonds
        int Nbonds;           // number of bonds

        int unitx;  // unit cell width and height (in terms of spins)
        int unity;  
        int x;      // lattice width and height (in terms of unit cells)
        int y;      

};



/******************************************************************************
 * Rectangular lattice class with periodic and open boundary conditions
 *****************************************************************************/
class Rectangle: public Lattice
{
    public:
        Rectangle(int _x, int _y, bool _OBC, long _seed);
};



/******************************************************************************
 * Chimera lattice class as implemented by Dwave. 
 * It assumes that the width of its unit cell is always kept to 2 spins
 *****************************************************************************/
class Chimera: public Lattice
{
    public:
        Chimera(int _x, int _y, long _seed,  int _unity);
};

#endif
