#include "bond.h"


/*****************************************************************************
 * Add a new bond 
 *****************************************************************************/
void Bonds::setBond(int _siteA, int _siteB)
/*
    siteA, siteB - coordinates of spins belonging to the newly defined bond  
*/
{
    tBond tbond; 
    bonds.push_back(tbond.set(_siteA,_siteB)); 
    Nbonds++;
}

/*****************************************************************************
 * Add a new bond with arbitrary strength value
 *****************************************************************************/
void Bonds::setBond(int _siteA, int _siteB, float _strength)
/*
    siteA, siteB - coordinates of spins belonging to the newly defined bond  
*/
{
    tBond tbond; 
    bonds.push_back(tbond.set(_siteA, _siteB, _strength)); 
    Nbonds++;
}


/*****************************************************************************
 * Return coordinates of spins associated with a bond 
 *****************************************************************************/
pair<int,int> Bonds::getSites(int index)
/*
    index - index of the bond 
*/
{
    return bonds[index].get(); 
}


/*****************************************************************************
 * Return bond's strength 
 *****************************************************************************/
float Bonds::getStrength(int index)
/*
    index - index of the bond 
*/
{
    return bonds[index].getStrength(); 
}



/*****************************************************************************
 * Print out the list of bonds
 *****************************************************************************/
void Bonds::print()
{
    
    cout << "   Bonds: " ;
    for (int i=0; i!=Nbonds; i++){
         bonds[i].print();
         cout << " ";
    }
    cout << endl;
}


//*****************************************************************************
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//*****************************************************************************



/******************************************************************************
 * Rectangle constructor - builds a rectangular lattice with PBC or OBC with
 * bonds only between nearest neighbours. 
 *****************************************************************************/
Rectangle::Rectangle(int _x, int _y, bool _OBC, float _J):
    Bonds() // Initialize the parent lattice class
/*
    x   - lattice width. 
    y   - lattice height. Caution: if a one-dimensional chain is required, 
                                   set y, and not x, to 1!
   _OBC - controls boundary conditions. Set to false for PBC and to true for OBC. 
*/
{
    // Set up the lattice type name and initialize its geometry
    name  = "rectangle";
    x = _x;
    y = _y;
    
    // Define assisting variables used in the loop below
    int siteA;  // coordinates of two spins
    int siteB;  // belonging to a bond 
    
    // One dimensional systems require a special treatement
    if (y == 1){
        for (int i=0; i!=x-int(_OBC); i++){
            
            // set up a signle bond per site
            siteA = i;
            siteB = (i+1)%x; // integer arithmetic is required 
                             // to properly treat boundary bonds
            setBond(siteA, siteB, _J);
        }
    }
    else{
        for (int i=0; i!=_y-int(_OBC); i++){
            for (int j=0; j!=_x-int(_OBC); j++){
                
                // set up a horizontal bond
                siteA = i*x+j;
                siteB = i*x+(j+1)%x; // integer arithmetic is required 
                                     // to properly treat boundary bonds
                setBond(siteA, siteB, _J);
                
                // set up a vertical bond
                siteA = i*x+j;
                siteB = ((i+1)%y)*x + j; // integer arithmetic is required 
                                         // to properly treat boundary bonds
                setBond(siteA, siteB, _J);
            }
        }
    }
}

//*****************************************************************************
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//*****************************************************************************



/******************************************************************************
 * Chimera constructor - builds a Dwave Chimera graph type lattice.
 * In doing so, it implicitely defines the lattice's coordinate system.
 * The 0'th spin is located in the upper left corner. The next spin is located
 * directly underneath it, etc. Once the count reaches the bottom left spin in the
 * current unit cell, it switches to the the upper-most spin located in the next
 * column of the same unit cell. Once all spins in the unit cell are enumerated,
 * the count continues with the unit cell on the right to the initial one.
 * In this way, all unit cells in the first row are enumerated. Then, it continues
 * to label sites in the second row in the same way.

 * Here is an example of chimera graph labelling convention for 2x2 grid with
 * 2x2 sites within each unit cell:

 *  0  2    4  6  
 *  1  3    5  7
 * 
 *  8 10   12 14 
 *  9 11   13 15
 *****************************************************************************/
Chimera::Chimera(int _x, int _y, int _unity, float _J):
    Bonds() 
{
    // Name the lattice type and initialize its geometry
    name  = "chimera";
    x     = _x;
    y     = _y;
    unitx = 2;
    unity = _unity;

    // Define assisting variables used in the loop below
    int LTspin = 0;  // coordinates of unit cell's left, top-most site
    int siteA;  // coordinates of two spins
    int siteB;  // belonging to a bond 

    // Construct all bonds by parsing through each unit cell one-by-one
    for (int i=0; i!=y; i++){
        for (int j=0; j!=x; j++){
            
            // Determine lattice coordinates of upper,
            // left-most  spin in the current unit cell
            LTspin = (2*unity)*(x*i + j); 
            
            // Add internal bonds within the current unit cell
            for (int k=0; k!=unity; k++){
                siteA = LTspin + k;
                for (int l=0; l!=unity; l++){
                    siteB = LTspin + unity + l;
                    setBond(siteA, siteB, _J);
                }
            }
            // Add vertical cross-cell bonds from the current unit cell
            if (i != y-1){
                for (int k=0; k!=unity; k++){
                    siteA = LTspin + k;
                    siteB = LTspin + k + 2*unity*x; 
                    setBond(siteA, siteB, _J);
                }
            }
            // Add horizontal cross-cell bonds from the current unit cell
            if  (j != x-1){
                for (int k=0; k!= unity; k++){
                    siteA = LTspin + unity + k;
                    siteB = LTspin + unity + k + 2*unity;
                    setBond(siteA, siteB, _J);
                }
            }

        }
    }
}

