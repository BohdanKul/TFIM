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

        vector<long>*  getSxs(){return &aSxs;}
        vector<float>* getSzs(){return &aSzs;}
        vector<float>* getSzSzs(){return &aSzSzs;}

        vector<long>*  getSxs0s(){return &aSxs_0slice;}
        vector<long>*  getSzs0s(){return &aSzs_0slice;}
        vector<long>*  getSzSzs0s(){return &aSzSzs_0slice;}
        
        long&          getCount(){return counter;}

        string Iheader;  // header for bonds 
        string Zheader;  // header for z-field
        string Xheader;  // header for x-field
        
        string Iheader0s;  // header for bonds 
        string Zheader0s;  // header for z-field
        string Xheader0s;  // header for x-field
    private:
        long DotProduct(vector<long>& v1, vector<long>& v2);
        long counter;
        int  Nspins;
        int  Nbonds;
        vector<long>  aSxs;
        vector<float> aSzs;
        vector<float> aSzSzs;
        vector<long>  aSxs_0slice;
        vector<long>  aSzs_0slice;
        vector<long>  aSzSzs_0slice;
        Spins* spins; 
        Bonds* bonds;
};

#endif

