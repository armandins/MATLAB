#include <mex.h>  /* fe_matrix_mex.c */
void lists(double n4e[], double c4n[],
        double Vol_T[], double Grads_T[],
        int nE, int nC, int d,
        double I[], double J[], double X[]){
    int j, m, n, r, idx1, idx2, ctr;
    double val;
    ctr = 0;
    for (j=0; j<nE; j++){
        for (m=0; m<d+1; m++){       
            for (n=0; n<d+1; n++){
                I[ctr] = n4e[j+m*nE]; J[ctr] = n4e[j+n*nE];
                val = 0.0;
                for (r=0; r<d; r++){
                   idx1 = j*(d+1)+m+r*(d+1)*nE;
                   idx2 = j*(d+1)+n+r*(d+1)*nE;
                   X[ctr] += Vol_T[j]*Grads_T[idx1]*Grads_T[idx2];
                }
                ctr += 1;
            }
        }
    }
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		 const mxArray *prhs[]){
    double *n4e, *c4n, *Vol_T, *Grads_T;
    int nE, nC, d;
    double *I, *J, *X;
    if (nrhs != 4)
        mexErrMsgTxt("4 input arguments required");   
    nC      = mxGetM(prhs[0]);
    d       = mxGetN(prhs[0]);
    nE      = mxGetM(prhs[1]);
    c4n     = mxGetPr(prhs[0]);
    n4e     = mxGetPr(prhs[1]);
    Vol_T   = mxGetPr(prhs[2]);
    Grads_T = mxGetPr(prhs[3]);
    if (nlhs != 3)
        mexErrMsgTxt("3 output arguments required");   
    plhs[0] = mxCreateDoubleMatrix(nE*(d+1)*(d+1),1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nE*(d+1)*(d+1),1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nE*(d+1)*(d+1),1,mxREAL);   
    I = mxGetPr(plhs[0]);
    J = mxGetPr(plhs[1]);
    X = mxGetPr(plhs[2]);
    lists(n4e,c4n,Vol_T,Grads_T,nE,nC,d,I,J,X);
}

