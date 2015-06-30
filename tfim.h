#ifndef TFIM_H
#define TFIM_H

#include <stack>
#include <list>
#include "communicator.h"
#include "hamiltonian.h"

using namespace std;


class tOperator
{
    public:
        int index; // bond or site index (depending on the operator type)
        int type;  // operator type; can take the following values:
                   // -2 null operator 
                   // -1 off-diagonal, one site
                   // 0 diagonal (identity), one site
                   // 1 diagonal, two sites
        tOperator(){type = -2; index = -1;};
        tOperator(int _type, int _index): type(_type), index(_index) {};
        tOperator set(int _type, int _index){type=_type; index=_index; return *this;};
};

class TFIM: public RandomBase
{
    public:
        TFIM(Spins* const _spins, Bonds* const _bonds, vector<float>* _xfield, long _seed, float _beta, long _binSize);
        int  DiagonalMove();       // diagonal operators update
        int  OffDiagonalMove();    // cluster update
        void MapStateBack();       // update of spins and operators list state
                                   // after an off-diagonal update 
        void ConstructLinks();     // construction of links required for the 
                                   // off-diagonal update
        void AdjustM();            // adjust the operator list length 
        void Measure();            // accumulate estimator measurements
        void Measure(int sampleInd, double* aEnergy, double* aMagnetization); // accumulate estimator measurements and populate output arrays

    private:
        void resetMeas();         // reset measurement variables to default values
        void computeDiagProb();   // compute diagonal operators insertation probabilities
        int VertexType(int otype, int index, Spins & ap); // return vertex type based on operator
                                              // type and its index in the list 
        int VtxToOperator(int vtxType);       // return operator type of a vertex type  
        int VtxFlip(int vtxType, int leg);    // return new vertex type after flipping
        
        void printIntVector(vector<int>* vect, string name );
        void printOperators();
       
        // Assististing datastructures in the diagonal update
        vector<float> diagProb;  // cumulative probability of a diagonal element 
        float bDiagProb;         // cumulative insertation probability  
                                 // of all bonds diagonal operators


        // Accumulate measurements for averaging 
        void Accumulate(long from,  long& to);
        void Accumulate(float from, float& to);
        void Accumulate(vector<float>& from, vector<float>& to);

        // Variables and data structures
        long  M;      // operators list length  
        long  n;      // number of non-identity operators
        float beta;   // inverse temperature
        int   Nbonds; // number of bonds
        int   Nsites; // number of sites

        vector<tOperator> lOper;     // operators list
        vector<int>       lVtx;      // vertices list
        vector<int>       lFirst;    // first and last leg coordinates
        vector<int>       lLast;     // for a particular spin site
        vector<int>       lLinks;    // linked vertex list
        vector<int>       lMarked;   // marked legs in a cluster update


        int LegSpin[8][4];  // vertex types to spin states map
        Hamiltonian ham;    // class incapsulating spins state and interactions
        //Spins   ap;       // propagated in imagenary time spins object


        // Accumulative variables for measurements
        vector<float> aMperSite; // magnetization per site
        vector<long> tMperSite;  // magnetization per site (assisting variable)
        long         an;         // operators list length
        float        aTM;        // total average magnetization
        long         nMeas;      // number of performed measurements

        // Flags to measure a particular estimator
        bool bMperSite; 
        bool bM;
        
        Communicator communicator;  // output files management
        long ID;
        long binSize;
        bool debug;
};
#endif
