#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <fstream>
#include <vector>
#include <unordered_map>
#include "lattice.h"

using namespace std;

class Communicator
{
    public:
        Communicator(Lattice * const lattice, float _beta, float _h, long _p);
        fstream* stream(string _fileName); 
        long     getId(){return id;};

        string dataName;
        string outDir;
        vector <string>  types;
        vector <fstream> files;   

    private:
        long id;
        long p;
        void GenerateId();
        unordered_map <string,fstream*> mFStreams;
}; 
#endif
