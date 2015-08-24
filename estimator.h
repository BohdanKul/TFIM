#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include "spin.h"
#include "bond.h"
#include <boost/format.hpp>

using namespace std;

class Estimator
{
    public:
        Estimator(Spins* const _spins, Bonds* const _bonds);
        void Accumulate(vector<vector<long>>& worldLines, long n);
        void Reset();    
        void print();
        string Iheader;  // header for bonds 
        string Zheader;  // header for z-field
        string Xheader;  // header for x-field
    private:
        long DotProduct(vector<long>& v1, vector<long>& v2);
        long counter;
        int  Nspins;
        int  Nbonds;
        vector<long>  aSxs;
        vector<float> aSzs;
        vector<float> aSzSzs;
        Spins* spins; 
        Bonds* bonds;
};

#endif

