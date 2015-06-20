#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <fstream>
#include <unordered_map>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include "bond.h"

using namespace std;

class Communicator
{
    public:
        Communicator(Bonds * const _bonds, float _beta, long _p);
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
