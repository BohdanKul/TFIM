#include <iostream>
#include <fstream>
#include "communicator.h"
#include <time.h>

using namespace std;
//using boost::lexical_cast;

Communicator::Communicator(Bonds* const _bonds, float _beta, long _p)
{
    p = _p;
     
    // Determine the filename 
    if (_bonds->getName()=="chimera")
       dataName = boost::str(boost::format("01-%03d-%03d-b%06.3f-x%02d-y%02d") 
                      %_bonds->getWidth() %_bonds->getHeight() %_beta%_bonds->getUnitWidth()%_bonds->getUnitHeight());
    
    if (_bonds->getName()=="rectangle")
       dataName = boost::str(boost::format("01-%03d-%03d-b%06.3f") 
                      %_bonds->getWidth() %_bonds->getHeight() %_beta);
    

    //Add the seed value to the filename
    dataName +=boost::str(boost::format("-p%05d") %p);

    types  = vector<string> {"estimator"};
    outDir = "OUTPUT"; 
    
    //Generate id
    GenerateId();
    string  fileName;
    for (vector<string>::iterator type=types.begin(); type!=types.end(); type++){
        fileName = boost::str(boost::format("%s/%s-%s-%09d.dat") %outDir %*type %dataName %id);
        mFStreams[*type] = new fstream(fileName,ios_base::out|ios_base::app);
        if (!*mFStreams[*type]){
           cerr << "Unable to process file: " << fileName << endl;
           exit(EXIT_FAILURE); 
        }
    }
    
}

void Communicator::GenerateId()
{
    time_t seconds = long(time(NULL) - 39*365*24*60*60);
    id = long(seconds + p);
    string fName;
    fName = boost::str(boost::format("OUTPUT/estimator-%s-%09d.dat") % dataName %id);
    boost::filesystem::path logPath(fName);

    while(boost::filesystem::exists(logPath)) {
        id += 1;
        fName = boost::str(boost::format("OUTPUT/estimator-%s-%09d.dat") % dataName %id);
        logPath = boost::filesystem::path(fName);

    }

}

fstream* Communicator::stream(string _fileName)
{
    return mFStreams[_fileName];

}


