
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
