#ifndef TFIM_H
#define TFIM_H

#include <stack>
#include <list>
#include "communicator.h"
#include "spin.h"

using namespace std;


class tOperator
{
    public:
        int type;  // operator type; can take the following values:
                   // 2  null operator 
                   // 1  diagonal
                   // 0  off-diagonal, 0'th leg flipped
                   // -1 off-diagonal, 1'th leg flipped 
        int index; // bond index
        tOperator(){type = 2; index = -1;};
        tOperator(int _type, int _index): type(_type), index(_index) {};
        tOperator set(int _type){type=_type; return *this;};
        tOperator set(int _type, int _index){type=_type; index=_index; return *this;};
};

class TFIM: public RandomBase
{
    public:
        TFIM(Spins* const _spins, Bonds* const _bonds, vector<float>* _xfield, vector<float>* _zfield, long _seed, float _beta, long _binSize);
        int  DiagonalMove();       // diagonal operators update
        int  OffDiagonalMove();    // cluster update
        void MapStateBack();       // update of spins and operators list state
                                   // after an off-diagonal update 
        void ConstructLinks();     // construction of links required for the 
                                   // off-diagonal update
        void AdjustM();            // adjust the operator list length
        bool AdjustLoopsN();       // adjust the number of built loops
        void Measure();            // accumulate estimator measurements
        void Measure(int sampleInd, double* aEnergy, double* aMagnetization); // accumulate estimator measurements and populate output arrays

    private:
        Spins& spins;
        Bonds& bonds;
        
        void resetMeas();         // reset measurement variables to default values
        void computeDiagProb();       // compute diagonal operators insertation probabilities
        
        void printIntVector(vector<int>* vect, string name );
        void printOperators();
       
        // Assististing datastructures in the diagonal update
        float  tdiagOffset;       // constant added to make all diagonal weigths non-negative
        vector<float> diagProb;  // cumulative probability of a diagonal element 

        vector<array<float,  4>> dWeights;     // weigths of diagonal vertices at each bond
 
        // Assististing datastructures in the -offdiagonal update
        vector<array<
                    array<
                         array<float,8>, 
                     16*8>, 
                  16>
              > VProb; // probability of exiting at a vertex leg
        // Accumulate measurements for averaging 
        void Accumulate(long from,  long& to);
        void Accumulate(float from, float& to);
        void Accumulate(vector<float>& from, vector<float>& to);

        // Variables and data structures
        long  M;      // operators list length  
        long  n;      // number of non-identity operators
        float beta;   // inverse temperature
        int   Nbonds; // number of bonds
        int   Nspins; // number of spins

        vector<tOperator> lOper;     // operators list
        vector<int>       lVtx;      // vertices list
        vector<int>       lFirst;    // first and last leg coordinates
        vector<int>       lLast;     // for a particular spin site
        vector<int>       lLinks;    // linked vertex list
        vector<int>       lBond;     // bonds associated with vertices
        vector<array<bool,5>> lClamped;

        int LegSpin[8][4];  // vertex types to spin states map
        //Spins   ap;       // propagated in imagenary time spins object


        // Accumulative variables for measurements
        float tDiagOffset;
        vector<float> aMperSite; // magnetization per site
        vector<long> tMperSite;  // magnetization per site (assisting variable)
        long         an;         // operators list length
        float        aTM;        // total average magnetization
        long         nMeas;      // number of performed measurements
        long  nFlippedLegs;      // number of flipped legs in the off-diagonal update
        long  nBounces;          // number of bounces in the off-diagonal update
        long  nLoops;
        float totalWeight;

        // Flags to measure a particular estimator
        bool bMperSite; 
        bool bM;
        
        Communicator communicator;  // output files management
        long ID;
        long binSize;
        bool debug;
};
#endif
