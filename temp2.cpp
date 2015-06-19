#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void fReadBond(fstream* bfile, int& siteA, int& siteB, float& strength)
{
    string         sbuf;
    istringstream ssbuf;
    getline(*bfile, sbuf);
    ssbuf.str(sbuf);
    ssbuf >> siteA;
    ssbuf >> siteB;
    ssbuf >> strength;
}

main()
{

    fstream Hamil ("hamil.dat", ios_base::in);
    string         sbuf;
    istringstream  ssbuf;
    getline(Hamil, sbuf);
    getline(Hamil, sbuf);
    ssbuf.str(sbuf);
    int nSz;
    int nSx;
    int nSzSz;
    ssbuf >> nSz;
    ssbuf >> nSx;
    ssbuf >> nSzSz;
    cout << " # Sz " << nSz << " # Sx " << nSx << " # SzSz " << nSzSz << endl; 
    int   siteA;
    int   siteB;
    float strength;
    cout << "---Sz: " << endl;
    for (int i=0; i!=nSz; i++){
        fReadBond(&Hamil, siteA, siteB, strength);
        cout << siteA << " " << siteB << " " << strength << endl;
    }
    cout << "---Sx: " << endl;
    for (int i=0; i!=nSx; i++){
        fReadBond(&Hamil, siteA, siteB, strength);
        cout << siteA << " " << siteB << " " << strength << endl;
    }
    cout << "---SzSz: " << endl;
    for (int i=0; i!=nSzSz; i++){
        fReadBond(&Hamil, siteA, siteB, strength);
        cout << siteA << " " << siteB << " " << strength << endl;
    }
}


