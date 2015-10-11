#include <boost/python.hpp>
#include <iostream>
#include "communicator.h"
#include "helper.h"
#include "vertex.h"
#include "directedloop.h"
#include "spin.h"
#include "bond.h"
#include "randombase.h"
#include "tfim.h"
#include "estimator.h"

namespace py = boost::python;
using namespace std;

// ------------------------------------------------------------------------ 
// Wrapper class around TFIM  
// ------------------------------------------------------------------------ 
class tWrapTFIM{
      public:
            tWrapTFIM(py::list dspins, py::list lbonds, py::list Z, py::list X, py::list ZZ, int N, float _beta, long seed);
           ~tWrapTFIM(); 
            int   DMove();
            void ODMove();
            void Adjust();
            void Measure(py::list Zs,  py::list Xs,  py::list ZZs); 
            void Measure(py::list E,   py::list n, 
                         py::list Zs,  py::list Xs,  py::list ZZs, 
                         py::list Z0s, py::list X0s, py::list ZZ0s);
            Spins* spins;
            Bonds* bonds;
            TFIM*  tfim;

            py::list Xfield;
            float beta;
            int NSpins;
};

// ------------------------------------------------------------------------ 
// Constructor
// ------------------------------------------------------------------------ 
tWrapTFIM::tWrapTFIM(py::list dspins, py::list lbonds, py::list Z, py::list X, py::list ZZ, int N, float _beta, long seed){

    NSpins = N;
    beta   = _beta;
    Xfield = X;

    // Initialize spins object
    if (len(dspins)>0){
        vector<int> vspins;
        for (int i=0; i<len(dspins); i++)
            vspins.push_back(py::extract<int>(dspins[i]));
            
        spins  = new Spins(N, seed, vspins);
    }
    else{
        spins  = new Spins(N, seed);
    }  
    
    // Initiate bonds object 
    bonds = new Bonds();
    int siteA; int siteB;
    float strength;
    cout << "Bonds: " << endl; 
    for (int i = 0; i < len(lbonds) ; i++){
        py::tuple b = py::extract<py::tuple>(lbonds[i]);
        siteA    = py::extract<int>(b[0]);
        siteB    = py::extract<int>(b[1]);
        strength = py::extract<float>(ZZ[i]);
        cout << "   " << siteA << " " << siteB << " " << strength << endl;
        bonds->setBond(siteA, siteB, strength);
    }

    vector<float> Zs;
    cout << "Zs: " ;
    for (int i = 0; i < len(Z) ; i++){
        Zs.push_back(py::extract<float>(Z[i]));
        cout << Zs[i] << " ";
    }
    cout << endl;
    
    vector<float> Xs;
    cout << "Xs: ";
    for (int i = 0; i < len(X) ; i++){
        Xs.push_back(py::extract<float>(X[i]));
        cout << Xs[i] << " ";
    }
    cout << endl;
    
    vector<float> ZZs;
    for (int i = 0; i < len(ZZ) ; i++){
        ZZs.push_back(py::extract<float>(ZZ[i]));
    }

    // Initiate MC object 
    int binSize = 100;
    tfim   = new TFIM(spins, bonds, &Xs, &Zs, seed, beta, binSize);

}

// --------------------------------------------------------------------------
// Measure averages of interest. Accumulating variables are reset in the end.
// -------------------------------------------------------------------------- 
void tWrapTFIM::Measure(py::list Zs, py::list Xs, py::list ZZs){
    
    Estimator* MCest = tfim->getEstimator();
    float nMeas = 1.0*MCest->getCount();

    int i=0;
    for (auto &Sx: *(MCest->getSxs())){
         Xs.append(-1.0*Sx/(nMeas*Xfield[i]*beta));
         i++;
    } 

    for (auto &Sz: *(MCest->getSzs())){
         Zs.append(1.0*Sz/nMeas);
    } 

    for (auto &SzSz: *(MCest->getSzSzs())){
         ZZs.append(1.0*SzSz/nMeas);
    } 

    tfim->resetMeas();
}


// --------------------------------------------------------------------------
// Measure averages of interest. Accumulating variables are reset in the end.
// -------------------------------------------------------------------------- 
void tWrapTFIM::Measure(py::list E, py::list n, py::list Zs, py::list Xs, py::list ZZs, py::list Z0s, py::list X0s, py::list ZZ0s){
    
    Estimator* MCest = tfim->getEstimator();
    float nMeas = 1.0*MCest->getCount();

    n.append(MCest->getAn()/nMeas);
    E.append(-1.0*MCest->getAn()/(nMeas*NSpins)/beta+tfim->getTDiagOffset()/(1.0*NSpins));

    int i=0;
    for (auto &Sx: *(MCest->getSxs())){
         Xs.append(-1.0*Sx/(nMeas*Xfield[i]*beta));
         i++;
    } 

    for (auto &Sz: *(MCest->getSzs())){
         Zs.append(1.0*Sz/nMeas);
    } 

    for (auto &SzSz: *(MCest->getSzSzs())){
         ZZs.append(1.0*SzSz/nMeas);
    } 

    i=0;
    for (auto &Sx: *(MCest->getSx0s())){
         X0s.append(-1.0*Sx/(nMeas*Xfield[i]*beta));
         i++;
    } 

    for (auto &Sz: *(MCest->getSz0s())){
         Z0s.append(1.0*Sz/nMeas);
    } 

    for (auto &SzSz: *(MCest->getSzSz0s())){
         ZZ0s.append(1.0*SzSz/nMeas);
    } 

    tfim->resetMeas();
}

// ------------------------------------------------------------------------ 
// Destructor 
// ------------------------------------------------------------------------ 
tWrapTFIM::~tWrapTFIM(){
    delete spins;
    delete bonds; 
    delete tfim;
} 

// ------------------------------------------------------------------------ 
// Diagonal move 
// ------------------------------------------------------------------------ 
int tWrapTFIM::DMove(){
    return tfim->DiagonalMove();
}

// ------------------------------------------------------------------------ 
// Parameters adjust  
// ------------------------------------------------------------------------ 
void tWrapTFIM::Adjust(){
    tfim->AdjustM();
}

// ------------------------------------------------------------------------ 
// Off-diagonal move
// ------------------------------------------------------------------------ 
void tWrapTFIM::ODMove(){
    tfim->ConstructLinks();
    tfim->OffDiagonalMove();
    tfim->MapStateBack();
}


// ============================================================================ 
// ============================================================================ 
BOOST_PYTHON_MODULE(tfim){
    py::class_<tWrapTFIM>("TFIM", py::init<py::list, py::list, py::list, py::list, py::list, int, float, long>())
        .def("DMove",   &tWrapTFIM::DMove)
        .def("ODMove",  &tWrapTFIM::ODMove)
        .def("Adjust",  &tWrapTFIM::Adjust)
        .def("Measure", &tWrapTFIM::Measure)
    ;
}
