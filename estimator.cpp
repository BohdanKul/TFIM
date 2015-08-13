#include "estimator.h"
#include "spin.cpp"
#include "bond.cpp"

/*****************************************************************************
 * Compute the dot product of a vector of +/- 1's with a shifted copy of itself.
 * To accelerate the execution, the vector's form is assumed to be a sequence
 * of integer values, each one representing the number of uninterrupted occurences
 * of +1's (-1's).
 *****************************************************************************/
long ShiftedDotProduct(vector<long>& v, long m){
    long tsum  = 0;   // the value of the dot product

    int sign = v[0]/abs(v[0]);  // sign of the current contribution 
    
    // Determine the index of the first block and its size 
    // in the shifted version of the original vector
    long shift = -m; // the block size
    long ishift = 0; // its index   
    shift = shift + v[0];
    while (!(shift>0)){
        ishift++;
        shift += abs(v[ishift]);
        sign = -1*sign;
    }
    sign = sign * v[0]/abs(v[0]);


    long blockL1 = abs(v[0]); // the size of the current block in the vector
    long blockL2 = shift;     // --------------------in the  shifted version  
    long bindex1 = 0;         // the index of the current block  in the vector
    long bindex2 = ishift;    // --------------------in the  shifted version  

    // Perfom the summation in the dot product.
    // First, exhaust entries in the shifted vector
    do{
        if (blockL1 < blockL2) {
            tsum += sign*blockL1;
            sign *= -1;
            
            blockL2 -= blockL1;
            
            bindex1 += 1;
            blockL1  = abs(v[bindex1]); 
        }
        else if (blockL1 > blockL2){
            tsum += sign*blockL2;
            sign *= -1;
            
            blockL1 -= blockL2;
            
            bindex2 += 1;
            if (bindex2!=v.size()) blockL2  = abs(v[bindex2]);
            else                   sign *= -1; 
        }
        else{
            tsum += sign*blockL1;
            bindex1 += 1;
            bindex2 += 1;
            blockL1  = abs(v[bindex1]); 
            if (bindex2!=v.size()) blockL2  = abs(v[bindex2]); 
            else                   sign *= -1; 
        }

    } while (bindex2 < v.size());

    // Switch to the start of the shifted vector
    blockL2 = v[0];
    bindex2 = 0;

    // Go through remaining entries of the original vector
    while (bindex1 < v.size()){
        if (blockL1 < blockL2) {
            tsum += sign*blockL1;
            sign *= -1;
            
            blockL2 -= blockL1;
            
            bindex1 += 1;
            blockL1  = abs(v[bindex1]); 
        }
        else if (blockL1 > blockL2){
            tsum += sign*blockL2;
            sign *= -1;
            
            blockL1 -= blockL2;
            
            bindex2 += 1;
            blockL2  = abs(v[bindex2]); 
        }
        else{
            tsum += sign*blockL1;
            bindex1 += 1;
            bindex2 += 1;
            if (bindex1!=v.size()) blockL1  = abs(v[bindex1]); 
            blockL2  = abs(v[bindex2]); 
        }
    }
}


/*****************************************************************************
* Estimator constructor 
*****************************************************************************/
Estimator::Estimator(Spins* const _spins, Bonds* const _bonds){
    counter = 0;
    spins   = _spins;
    bonds   = _bonds;
    Nspins = spins->getSize();
    Nbonds = bonds->getBondsN();
    
    aSzs.resize(Nspins,0);
    aSxs.resize(Nspins,0);
    aSzSzs.resize(Nbonds,0);

    Iheader = "";  // header for interactions
    Zheader = "";  // header for z-field
    Xheader = "";  // header for x-field

    // create header for the bonds
    int is0; int is1;
    for (auto ibond=0; ibond!=Nbonds; ibond++){
        is0 = bonds->getBond(ibond)->getSiteA();
        is1 = bonds->getBond(ibond)->getSiteB();
        if (ibond==0) Iheader += boost::str(boost::format("#(%4d,%4d)    ")%is0%is1); 
        else          Iheader += boost::str(boost::format(" (%4d,%4d)    ")%is0%is1); 
    }

    // create headers for the fields
    for (auto i=0; i!=Nspins; i++){
        if (i==0) Zheader += boost::str(boost::format("#%4d           ")%i); 
        else      Zheader += boost::str(boost::format(" %4d           ")%i); 
    }
    
    Xheader = Zheader;
}

/*****************************************************************************
* Accumulate measurements
*****************************************************************************/
void Estimator::Accumulate(vector<vector<long>>& worldLines, long n){
    counter += 1;
    for (auto ispin=0; ispin!=Nspins; ispin++){
        aSxs[ispin] += worldLines[ispin].size(); 
    }

    vector<long>* wl;
    long spinT;
    for (auto ispin=0; ispin!=Nspins; ispin++){
        wl = &(worldLines[ispin]);
        if ((wl->size()==0) or (wl->back() != n))
            wl->push_back(n);
        spinT = wl->at(0);
        for (auto fSlice=1; fSlice!=wl->size(); fSlice++){
            spinT += (wl->at(fSlice) - wl->at(fSlice-1))*pow(-1,fSlice);      
        }
        aSzs[ispin] += (1.0*spinT*spins->getSpin(ispin))/((float) n+1.0);
    }

    int  s0; int  s1;
    int is0; int is1;
    for (auto ibond=0; ibond!=Nbonds; ibond++){
        is0 = bonds->getBond(ibond)->getSiteA();
        is1 = bonds->getBond(ibond)->getSiteB();
        s0  = spins->getSpin(is0);
        s1  = spins->getSpin(is1);
        aSzSzs[ibond] += (1.0*s0*s1*DotProduct(worldLines[is0], worldLines[is1]))/((float) n+1.0); 
    }
}


/*****************************************************************************
* Reset accumulating variables to zeros 
*****************************************************************************/
void Estimator::Reset(){    
    counter = 0;
    fill(aSxs.begin(),   aSxs.end(),   0);
    fill(aSzs.begin(),   aSzs.end(),   0);
    fill(aSzSzs.begin(), aSzSzs.end(), 0);
}


/*****************************************************************************
 * Compute the dot product of two vectors of +/- 1's. To accelerate the operation, 
 * the vector's form is assumed to be a sequence of integer values, each one 
 * representing the number of uninterrupted occurences of +1's (-1's).
 *****************************************************************************/
long Estimator::DotProduct(vector<long>& v1, vector<long>& v2){
    long blockL1 = v1[0]; // the size of the current block in the vector
    long blockL2 = v2[0]; 
    long bindex1 = 0;     // current block
    long bindex2 = 0;

    int sign = 1;   // alternating sign
    int tsum = 0;   // total sum
    // Perfom the summation in the dot product.
    do{
        if (blockL1 < blockL2) {
            tsum += sign*blockL1;
            sign *= -1;
            blockL2 -= blockL1;
            
            bindex1 += 1;
            if (bindex1!=v1.size()) blockL1  = v1[bindex1] - v1[bindex1-1]; 
        }
        else if (blockL1 > blockL2){
            tsum += sign*blockL2;
            sign *= -1;
            
            blockL1 -= blockL2;
            
            bindex2 += 1;
            if (bindex2!=v2.size()) blockL2  = v2[bindex2] - v2[bindex2-1];
        }
        else{
            tsum += sign*blockL1;
            bindex1 += 1;
            bindex2 += 1;
            if (bindex1!=v1.size()) blockL1  = v1[bindex1] - v1[bindex1-1]; 
            if (bindex2!=v2.size()) blockL2  = v2[bindex2] - v2[bindex2-1]; 
        }

    } while ((bindex1 < v1.size()) and (bindex2 < v2.size()));
    
    return tsum;
}

/*****************************************************************************
* Print out accumulated variables 
*****************************************************************************/
void Estimator::print(){
    cout << "Szs: " << endl;
    for (int i=0; i!=Nspins; i++)
        cout << aSzs[i] << " ";

    cout << endl << "Sxs: " << endl;
    for (int i=0; i!=Nspins; i++)
        cout << aSxs[i] << " ";

    cout << endl << "SzSzs: " << endl;
    for (int i=0; i!=Nbonds; i++)
        cout << aSzSzs[i] << " ";
    cout << endl;
}

