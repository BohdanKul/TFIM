#include "lattice.h"

/*****************************************************************************
 * Basis lattice constructor
 *****************************************************************************/
Spins::Spins(int _Nspins)
/*
    _Nspins sets the size of spins vectors
*/
{
    // Initialize varialbes
    Nspins = _Nspins;
    spins.resize(Nspins,1);
}

/*****************************************************************************
 * Basis lattice constructor - initializes spins in a random state
 *****************************************************************************/
Spins::Spins(int _Nspins, long _seed):
    Spins(_Nspins)
/*
    _Nspins - sets the size of spins vectors
    _seed   - initializes the random number generator necessary 
              to generate a random distribution of spins
*/
{
    RandomBase  rand(_seed);     // random numbers generator
    
    // Initialize spins randomly with values 1 and -1
    for (auto spin=spins.begin(); spin!=spins.end(); spin++) {                          
        *spin = pow(-1,rand.uRandInt()%2);                         
    }
}

/*****************************************************************************
 * Returns the spin state 
 *****************************************************************************/
int Spins::getSpin(int index)
/*
   index - index of the element to be flipped  
*/
{
    return spins[index];
}

/*****************************************************************************
 * Sets the spin state to a particular value
 *****************************************************************************/
void Spins::setSpin(int index, int val)
/*
   index - index of the element to be flipped  
   val   - the new value of the spin state
*/
{
    spins[index] = val;
}

/*****************************************************************************
 * Flips a spin
 *****************************************************************************/
void Spins::flipSiteSpin(int index)
/*
   index - index of the element to be flipped  
*/
{
    setSpin(index, -1*getSpin(index));
}

/*****************************************************************************
 * Print out the spins state 
 *****************************************************************************/
void Spins::printSpins()
{
    cout << "   Spins state: " << endl;
    for (int i=0; i!=Nspins; i++)
        cout << getSpin(i) << " ";
    cout << endl;
}

//*****************************************************************************
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//*****************************************************************************


/*****************************************************************************
 * Base lattice constructor.
 * Initiate spins in a random state and set up the basic geometry.
 *****************************************************************************/
Lattice::Lattice(int _x, int _y, int _unitx, int _unity, long _seed):
    Spins(_unitx*_unity*_x*_y, _seed)  // initialize the parental class
/*
    _seed - random seed generator
*/
{
    name = "general lattice";

    // Set up geometric parameters of the underlying lattice
    unitx = _unitx;
    unity = _unity;
    x     = _x;
    y     = _y;

    // Bonds are to be defined in daugther classes
    Nbonds = 0;
}

/*****************************************************************************
 * Add a new bond 
 *****************************************************************************/
void Lattice::setBond(int siteA, int siteB)
/*
    siteA, siteB - coordinates of spins belonging to the newly defined bond  
*/
{
    tBond tbond; 
    bonds.push_back(tbond.set(siteA,siteB)); 
    Nbonds++;
}

/*****************************************************************************
 * Return coordinates of spins associated with a bond 
 *****************************************************************************/
pair<int,int> Lattice::getSites(int index)
/*
    index - index of the bond 
*/
{
    return bonds[index].get(); 
}

/*****************************************************************************
 * Return the state of two spins associated with a bond
 *****************************************************************************/
pair<int,int> Lattice::getSpins(int index)
/*
    index - index of the bond 
*/
{
    pair<int,int> sites = getSites(index);
    return pair<int,int>  (getSpin(sites.first),getSpin(sites.second));  
}


/*****************************************************************************
 * Flip both spins associated with a bond 
 *****************************************************************************/
void Lattice::flipBondSpins(int index)
/*
    index - index of the bond 
*/
{
    pair<int,int> sites = getSites(index);
    flipSiteSpin(sites.first);
    flipSiteSpin(sites.second);
}

/*****************************************************************************
 * Print out the list of bonds
 *****************************************************************************/
void Lattice::printBonds()
{
    
    cout << "   Lattice bonds: " ;
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
Rectangle::Rectangle(int _x, int _y, bool _OBC = false, long _seed = 0):
    Lattice(_x,_y,1,1, _seed) // Initialize the parent lattice class
                              // with a unit cell containing only 1 spin.
/*
    x   - lattice width. 
    y   - lattice height. Caution: if a one-dimensional chain is required, 
                                   set y, and not x, to 1!
   _OBC - controls boundary conditions. Set to false for PBC and to true for OBC. 
*/
{
    // Set up the lattice type name
    name  = "rectangle";
    
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
            setBond(siteA,siteB);
        }
    }
    else{
        for (int i=0; i!=_y-int(_OBC); i++){
            for (int j=0; j!=_x-int(_OBC); j++){
                
                // set up a horizontal bond
                siteA = i*x+j;
                siteB = i*x+(j+1)%x; // integer arithmetic is required 
                                     // to properly treat boundary bonds
                setBond(siteA,siteB);
                
                // set up a vertical bond
                siteA = i*x+j;
                siteB = ((i+1)%y)*x + j; // integer arithmetic is required 
                                         // to properly treat boundary bonds
                setBond(siteA,siteB);
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
 * the count continues with the next unit cell located underneath the initial one.
 * The sequence of unit cells is the same snake-like sequence as the one used
 * within a unit cell.
 *****************************************************************************/
Chimera::Chimera(int _x, int _y, long _seed = 0, int _unity = 4):
    Lattice(_x, _y, 2, _unity, _seed) 
{
    // Name the lattice type
    name  = "chimera";

    // Define assisting variables used in the loop below
    int LTspin = 0;  // coordinates of unit cell's left, top-most site
    int siteA;  // coordinates of two spins
    int siteB;  // belonging to a bond 

    // Construct all bonds by parsing through each unit cell one-by-one
    for (int i=0; i!=y; i++){
        for (int j=0; j!=x; j++){
            
            // Determine lattice coordinates of upper,
            // left-most  spin in the current unit cell
            LTspin = (2*unity)*(i + y*j); 
            
            // Add internal bonds within the current unit cell
            for (int k=0; k!=unity; k++){
                siteA = LTspin + k;
                for (int l=0; l!=unity; l++){
                    siteB = LTspin + unity + l;
                    setBond(siteA,siteB);
                }
            }
            // Add vertical cross-cell bonds from the current unit cell
            if (i != y-1){
                for (int k=0; k!=unity; k++){
                    siteA = LTspin + k;
                    siteB = LTspin + k + 2*unity; 
                    setBond(siteA,siteB);
                }
            }
            // Add horizontal cross-cell bonds from the current unit cell
            if  (j != x-1){
                for (int k=0; k!= unity; k++){
                    siteA = LTspin + unity + k;
                    siteB = LTspin + unity + k + 2*unity*y;
                    setBond(siteA,siteB);
                }
            }

        }
    }
}

