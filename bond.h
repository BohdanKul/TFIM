#ifndef BOND_H
#define BOND_H

#include <string.h>
#include <vector>
#include <iostream>
using namespace std;




/******************************************************************************
 * Class holding coordinates of 2 sites belonging to a bond, 
 * as well as its strength value.
 *****************************************************************************/
class tBond
{
    public:

        tBond(): A(-1), B(-1), strength(-1) {};                  // initiate the bond to default values
        tBond(int x, int y): A(x), B(y), strength(-1) {};         // initiate bond's spins to user-defined indices  
        tBond(int x, int y, float s): A(x), B(y), strength(s) {}; // initiate bond's spins and bond's strength to 
                                                                  // user-defined values

        tBond set(int x, int y)            // reset bond's spins to new values
            { A = x; B = y; return *this;};
        tBond set(int x, int y, float s)   // reset bond's spins and value it to new values
            { A = x; B = y; strength = s; return *this;};
        pair<int,int> get()                // get the indices of bond's spins 
            {return pair<int,int>(A,B);};
        int getSiteA()                 // get the index bond's first spin 
            {return A};
        int getSiteB()                 // get the index bond's second spin 
            {return B};
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
 * Virtual class that serves as a starting point of implementation of bonds 
 * acting on two sites. This class only  implements the necessary functionality 
 * to deal with bonds without actually setting them up. In order to define 
 * physical bonds, one needs to inherit this class and define them in the 
 * constructors of the new class using setBond() method. For an example, 
 * see rectangle and chimera classes.
 *****************************************************************************/
class Bonds 
{
    public:
             
        Bonds(): name("general"), Nbonds(0), x(-1), y(-1), unitx(-1), unity(-1) {};

        void setBond(int _siteA, int _siteB);                  // set bond between two sites
        void setBond(int _siteA, int _siteB, float _strength);  // same but with a specific strength value
        
        void print();                        // print out all bonds 

        pair<int,int> getSites(int index);          // get sites associated with a bond
        float         getStrength(int index);       // get bond's strength
        tBond*        getBond(int index)
                      { return &(bonds[index]) }
        int           getBondsN(){ return Nbonds;}; // get the number of bonds
        string&       getName(){return name;};      // get the lattice name

        int getWidth(){      return x;};
        int getHeight(){     return y;};
        int getUnitWidth(){  return unitx;};
        int getUnitHeight(){ return unity;};
    protected:
        string name;          // lattice-type name
        
        vector<tBond> bonds;  // vector of bonds
        int Nbonds;           // number of bonds

        int x;      // width
        int y;      // height
        int unitx;  // unit cell width
        int unity;  // unit cell height
};



/******************************************************************************
 * Implements nearest neighbors bonds on a PBC or OBC rectangular lattice
 *****************************************************************************/
class Rectangle: public Bonds 
{
    public:
        Rectangle(int _x, int _y, bool _OBC, float _J);
};



/******************************************************************************
 * Implements chimera-type bonds as used on D-wave 2. 
 * It assumes that the width of its unit cell is always kept to 2 spins
 *****************************************************************************/
class Chimera: public Bonds
{
    public:
        Chimera(int _x, int _y, int _unity, float _J);
};


#endif
