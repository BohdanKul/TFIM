#include <mex.h>
#include "qmc.h"
//matlab usage:
//[indices] = unique_tol(M, cutoff);
void printSyntax() {
	mexPrintf("qmc(int* ind1, int* ind2, double* J, int nJs, int nSpins, double T, double Delta, int nSamples, int binSize, int seed)");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
*prhs[]) {
    if ((nrhs != 9) ) {
		printSyntax();
		mexErrMsgTxt("Invalid number of input arguments");
	}
    double* t1 = mxGetPr(prhs[0]);
    double* t2 = mxGetPr(prhs[1]);	
    double* J = mxGetPr(prhs[2]);
    int nJs = mxGetN(prhs[2]);	
    int ind1[nJs];
    int ind2[nJs];
    for(int i=0;i!=nJs;i++){
        ind1[i] = (int) t1[i] - 1;
        ind2[i] = (int) t2[i] - 1;
    }
    
	int nSpins = (int) mxGetScalar(prhs[3]);
    double T = mxGetScalar(prhs[4]);
    double Delta = mxGetScalar(prhs[5]);
    int nSamples = (int) mxGetScalar(prhs[6]);
    int binSize = (int) mxGetScalar(prhs[7]);
    int seed = (int) mxGetScalar(prhs[8]);	
    double aEnergy[nSamples]; 
    double aMagnetization[nSamples];
    qmc(ind1, ind2, J, nJs, nSpins, T, Delta, nSamples, binSize, seed, aEnergy, aMagnetization);	
    plhs[0] = mxCreateDoubleMatrix(1, nSamples, mxREAL);
	double * ae = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, nSamples, mxREAL);
	double * am = mxGetPr(plhs[1]);
	for (int i = 0; i != nSamples; i++) {
        *ae++ = aEnergy[i];
        *am++ = aMagnetization[i];
	}
}
