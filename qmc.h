#include <cstdlib>
#include <cmath>
#include <boost/format.hpp> 
#include <boost/program_options.hpp>
#include "tfim.h" 
void qmc(int* ind1, int* ind2, double* J, int nJs, int nSpins, double T, double Delta, int nSamples, int binSize, int seed, double* aEnergy, double* aMagnetization) {
     // ------------------------------------------------------------------------ 
    // Initialize interactions
    // ------------------------------------------------------------------------ 
    Bonds* bonds = new Bonds();
    vector<float> xfield(nSpins, Delta);
    for(int i=0;i<nJs; i++){       
        bonds->setBond(ind1[i], ind2[i], J[i]);
        cout << "   Bond (" << ind1[i] << "," << ind2[i] << ") = " << J[i] << endl;   
    }
     
    // ------------------------------------------------------------------------ 
    // Initialize spins 
    // ------------------------------------------------------------------------ 
    Spins * spins = new Spins(nSpins, seed);
    
    // ------------------------------------------------------------------------ 
    // Initialize Monte Carlo class 
    // ------------------------------------------------------------------------ 
    TFIM tfim(spins, bonds, &xfield, seed, 1/T, binSize);    

    // ------------------------------------------------------------------------ 
    // Run the main loop 
    // ------------------------------------------------------------------------ 
    cout << endl << "Measurement stage" << endl << endl;
    int sampleInd = 0;
    for (long i=0; i!=binSize*nSamples; i++){
        // Perform  a full MC sweep
        if((i+1) % binSize == 0)
            sampleInd = (i+1)/binSize;
            
        while (tfim.DiagonalMove()==1)               // diagonal update
              tfim.AdjustM();
              
        tfim.ConstructLinks();                       // linked list and vertex list construction 
        tfim.OffDiagonalMove();                      // cluster update
        tfim.MapStateBack();                         // mapping back the updated state 
        tfim.Measure(sampleInd, aEnergy, aMagnetization);
    }

    delete spins;
    delete bonds;   
}