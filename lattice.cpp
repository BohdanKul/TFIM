#include "lattice.h"

/*****************************************************************************
 * Basis lattice constructor
 *****************************************************************************/
Spins::Spins(int _Nsites)
{
    // Initialize varialbes
    Nsites = _Nsites;
    spins.resize(Nsites,0);
}

/*****************************************************************************
 * Basis lattice constructor - initializes spins in a random state
 *****************************************************************************/
Spins::Spins(int _Nsites, long _seed):
    Spins(_Nsites)
{
    RandomBase  rand(_seed);     // random numbers generator
    
    // Initialize spins randomly with values 1 and -1
    for (auto spin=spins.begin(); spin!=spins.end(); spin++) {                          
        *spin = pow(-1,rand.uRandInt()%2);                         
    }
}

/*****************************************************************************
 * Flips a spin indicated by its index 
 *****************************************************************************/
void Spins::flipSiteSpin(int index)
{
    setSpin(index, -1*getSpin(index));
}

/*****************************************************************************
 * Returns the spin state indicated by its index 
 *****************************************************************************/
int Spins::getSpin(int index)
{
    return spins[index];
}

/*****************************************************************************
 * Returns the spin state indicated by its index 
 *****************************************************************************/
void Spins::setSpin(int index, int val)
{
    spins[index] = val;
}

/*****************************************************************************
 * Print out the spins state 
 *****************************************************************************/
void Spins::printSpins()
{
    cout << "Spins state " << endl;
    for (int i=0; i!=Nsites; i++)
        cout << spins[i] << " ";
    cout << endl;
}

/*****************************************************************************
 * Base lattice constructor - initiates spins in a random state
 *****************************************************************************/
Lattice::Lattice(int _x, int _y, int _unitx, int _unity, long _seed):
    Spins(_unitx*_unity*_x*_y, _seed)
{
    // Initialize varialbes
    unitx = _unitx;
    unity = _unity;
    x     = _x;
    y     = _y;
    Nbonds = 0;
}

/*****************************************************************************
 * Adds a new bond 
 *****************************************************************************/
void Lattice::setBond(int siteA, int siteB)
{
    tBond tbond;                             
    bonds.push_back(tbond.set(siteA,siteB)); 
    Nbonds++;
}

/*****************************************************************************
 * Return lattice sites' indices associated with a bond
 *****************************************************************************/
pair<int,int> Lattice::getSites(int index)
{
    tBond bond = bonds[index];
    return pair<int,int> (bond.A, bond.B); 
}

/*****************************************************************************
 * Return the spins state associated with a bond
 *****************************************************************************/
pair<int,int> Lattice::getSpins(int index)
{
    tBond bond = bonds[index];
    return pair<int,int>  (getSpin(bond.A),getSpin(bond.B));  
}


/*****************************************************************************
 * Flip spins belonging to a bond indicated by its index 
 *****************************************************************************/
void Lattice::flipBondSpins(int index)
{
    tBond bond = bonds[index];
    flipSiteSpin(bond.A);
    flipSiteSpin(bond.B);
}

/*****************************************************************************
 * Print out the spins state 
 *****************************************************************************/
void Lattice::printBonds()
{
    cout << "   Bond: " ;
    for (int i=0; i!=Nbonds; i++)
        cout << "(" << bonds[i].A << "," << bonds[i].B << ") ";
    cout << endl;
}



/******************************************************************************
 * Rectangle constructor - builds a rectangular lattice with PBC or OBC
 *****************************************************************************/
Rectangle::Rectangle(int _x, int _y, bool _OBC = false, long _seed = 0):
    Lattice(_x,_y,1,1, _seed)
{
    // Define assisting variables used in the loop below
    int siteA;
    int siteB;
    name  = "rectangle";
    // One dimensional systems require a special treatement
    if (y == 1){
        for (int i=0; i!=x-int(_OBC); i++){
            siteA = i;
            siteB = (i+1)%x;
            setBond(siteA,siteB);
        }
    }
    else{
        for (int i=0; i!=_y-int(_OBC); i++){
            for (int j=0; j!=_x-int(_OBC); j++){
                siteA = i*x+j;
                siteB = i*x+(j+1)%x;
                setBond(siteA,siteB);
                
                siteA = i*x+j;
                siteB = ((i+1)%y)*x + j;
                setBond(siteA,siteB);
            }
        }
    }
}

/******************************************************************************
 * Chimera constructor - builds a Dwave Chimera graph type lattice
 * The sites' count first iterates through a unit cell starting in downwards
 * vertical direction. Then, it continues through each unit cell in the same way.
 *****************************************************************************/
Chimera::Chimera(int _x, int _y, long _seed = 0, int _unity = 4):
    Lattice(_x, _y, 2, _unity, _seed) 
{
    // Define assisting variables used in the loop below
    int LTspin = 0;  // coordinates of unit cell's left, top-most site
    int siteA;    
    int siteB;   
    name  = "chimera";

    // Construct all bonds by parsing through each unit cell one-by-one
    for (int i=0; i!=y; i++){
        for (int j=0; j!=x; j++){
            LTspin = (2*unity)*(i + y*j); 
            
            // Add internal bonds within a unit cell
            for (int k=0; k!=unity; k++){
                siteA = LTspin + k;
                for (int l=0; l!=unity; l++){
                    siteB = LTspin + unity + l;
                    setBond(siteA,siteB);
                }
            }
            // Add vertical cross-cell bonds from a unit cell
            if (i != y-1){
                for (int k=0; k!=unity; k++){
                    siteA = LTspin + k;
                    siteB = LTspin + k + 2*unity; 
                    setBond(siteA,siteB);
                }
            }
            // Add horizontal cross-cell bonds from a unit cell
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

