/*
[k1A k1B] = growAB(A, B);
[k2A k2B] = growAB(A+h2.*k1A, B+h2.*k1B);
[k3A k3B] = growAB(A+h2.*k2A, B+h2.*k2B);
[k4A k4B] = growAB(A+h.*k3A, B+h.*k3B);
A = A + h6 * (k1A + 2*k2A + 2*k3A + k4A);    
B = B + h6 * (k1B + 2*k2B + 2*k3B + k4B);

function [dA dB] = growAB(A, B)
    dA = m*B.*(1-A);
    dB = p*A - dA;
end
*/
   

/* $Revision: 1.10.6.4 $ */
#include <math.h>
#include "mex.h"

/* Input Arguments */
#define	A_IN	prhs[0]
#define	B_IN	prhs[1]
#define h_IN    prhs[2]
#define r_IN    prhs[3]
#define p_IN    prhs[4]

/* Output Arguments */
#define	A_OUT	plhs[0]
#define	B_OUT	plhs[1]

size_t numel;

static void growAB_RK4(double A[], double B[], double h, double r, double p, double oA[], double oB[]) {
    double dA[4], dB[4];
    double tA, tB;
    double h2 = h/2;
    double h6 = h/6;
    int i;
 
    for(i=0; i<numel; i++) {
        dA[0] = r*B[i] * (1-A[i]);
        dB[0] = p*A[i] - dA[0];
        
        tA = A[i]+h2*dA[0];
        tB = B[i]+h2*dB[0];        
        dA[1] = r*tB * (1-tA);
        dB[1] = p*tA - dA[1];
        
        tA = A[i]+h2*dA[1];
        tB = B[i]+h2*dB[1];
        dA[2] = r*tB * (1-tA);
        dB[2] = p*tA - dA[2];        
        
        tA = A[i]+h*dA[2];
        tB = B[i]+h*dB[2];
        dA[3] = r*tB * (1-tA);
        dB[3] = p*tA - dA[3];
        
        oA[i] +=  h6*(dA[0] + 2*dA[1] + 2*dA[2] + dA[3]);
        oB[i] +=  h6*(dB[0] + 2*dB[1] + 2*dB[2] + dB[3]);
    }    
    
    return;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) { 
    double *A, *B, h, r, p; 
    double *oA, *oB; 
    int nA, nB;
    
    /* Check for proper number of arguments */
    if (nrhs != 5) { 
        mexErrMsgTxt("Five input arguments required."); 
    } else if (nlhs != 2) {
        mexErrMsgTxt("Two outputs required."); 
    } 
    
    /* Check the dimensions of A, B */ 
    nA = mxGetNumberOfElements(A_IN);
    nB = mxGetNumberOfElements(B_IN);
    numel = nA;
    if (!mxIsDouble(A_IN) || mxIsComplex(A_IN) || !mxIsDouble(B_IN) || mxIsComplex(B_IN) || (nA != nB)) { 
        mexErrMsgTxt("A & B must be real matrices that are the same size."); 
    } 
    
    /* Create a matrix for the return argument */ 
    A_OUT = mxDuplicateArray(A_IN); 
    B_OUT = mxDuplicateArray(B_IN);
    
    /* Assign pointers to the various parameters */ 
    A = mxGetPr(A_IN);
    B = mxGetPr(B_IN);
    r = mxGetScalar(r_IN);
    p = mxGetScalar(p_IN);
    h = mxGetScalar(h_IN);
    oA = mxGetPr(A_OUT);
    oB = mxGetPr(B_OUT);
        
    /* Do the actual computations in a subroutine */
    growAB_RK4(A,B,h,r,p,oA,oB); 
    return;
}
