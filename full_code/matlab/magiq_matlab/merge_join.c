/*==========================================================
 * merge_join.c - join two sorted arrays with duplicates
 *
 * Merges two sorted arrays (A and B) and produces two index arrays 
 * (IA and IB), where A(IA) = B(IB).
 *
 * The calling syntax is:
 *
 *		[ IA  IB ] = merge_join(A, B)
 *
 * This is a MEX-file for MATLAB.
 * Fuad Jamour ;)
 * Algorithm from: [http://www.dcs.ed.ac.uk/home/tz/phd/thesis/node20.htm]
 *========================================================*/

#include "mex.h"
#include "string.h"

mwSize
expand_out
(
    mxArray  **IA_mx_ptr,
    mxArray  **IB_mx_ptr,
    mwSize   o_sz
)
{
    //mwSize new_sz = o_sz * 2;
    mwSize new_sz = o_sz * 20;
    
    /* create new arrays, copy data, then destroy old arrays */
    mxArray *IA_mx_ptr_new = mxCreateDoubleMatrix(new_sz, 1, mxREAL);
    mxArray *IB_mx_ptr_new = mxCreateDoubleMatrix(new_sz, 1, mxREAL);
    
    memcpy(mxGetPr(IA_mx_ptr_new), mxGetPr(*IA_mx_ptr), sizeof(double) * o_sz);
    memcpy(mxGetPr(IB_mx_ptr_new), mxGetPr(*IB_mx_ptr), sizeof(double) * o_sz);
    
    mxDestroyArray(*IA_mx_ptr);
    mxDestroyArray(*IB_mx_ptr);
    
    *IA_mx_ptr = IA_mx_ptr_new;
    *IB_mx_ptr = IB_mx_ptr_new;
    return new_sz;
}

void
shrink_out
(
    mxArray  **IA_mx_ptr,
    mxArray  **IB_mx_ptr,
    mwSize   o_sz
)
{
    double *ptrA    = mxGetPr(*IA_mx_ptr);
    double *ptrB    = mxGetPr(*IB_mx_ptr);
    void   *newptrA = mxRealloc(ptrA, sizeof(double) * o_sz);
    void   *newptrB = mxRealloc(ptrB, sizeof(double) * o_sz);
    mxSetPr(*IA_mx_ptr, newptrA);
    mxSetPr(*IB_mx_ptr, newptrB);
    mxSetM(*IA_mx_ptr, o_sz);
    mxSetM(*IB_mx_ptr, o_sz);
}

mwSize
out_set
(
    double **IAp,
    double **IBp,
    mxArray  **IA_mx_ptr,
    mxArray  **IB_mx_ptr,
    mwSize i,
    mwSize j,
    mwSize o,
    mwSize o_sz
)
{
    double *IA = *IAp;
    double *IB = *IBp;
    IA[o] = i+1; IB[o] = j+1;
    if(o == o_sz-1) {
        o_sz = expand_out(IA_mx_ptr, IB_mx_ptr, o_sz);
        *IAp = mxGetPr(*IA_mx_ptr);
        *IBp = mxGetPr(*IB_mx_ptr);
    }
    return o_sz;    
}

void
merge_join_double
(
    double   *A,
    double   *B,
    mxArray  **IA_mx_ptr,
    mxArray  **IB_mx_ptr,
    mwSize   A_size, 
    mwSize   B_size)
{
    mwSize i = 0;
    mwSize j = 0;
    mwSize o = 0;
    
    mwSize ii;
    mwSize jj;
    
    /* allocate output */
    mwSize o_sz = A_size + B_size;
    *IA_mx_ptr = mxCreateDoubleMatrix(o_sz, 1, mxREAL);
    *IB_mx_ptr = mxCreateDoubleMatrix(o_sz, 1, mxREAL);
    
    double *IA = mxGetPr(*IA_mx_ptr);
    double *IB = mxGetPr(*IB_mx_ptr);
    
    
    while(i != A_size && j != B_size) {
        if     (A[i] > B[j]) ++j;
        else if(A[i] < B[j]) ++i;
        else {//A[i] == B[j]
            o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, i, j, o, o_sz);
            ++o;
            jj = j+1;
            while(jj != B_size && A[i] == B[jj]) {
                o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, i, jj, o, o_sz);
                ++o;
                ++jj;
            }
            ii = i+1;
            while(ii != A_size && B[j] == A[ii]) {
                o_sz = out_set(&IA, &IB, IA_mx_ptr, IA_mx_ptr, ii, j, o, o_sz);
                ++o;
                ++ii;
            }
            ++i; ++j;
        }
    }
    
    /* resize the output to avoid trailing zeros */
    shrink_out(IA_mx_ptr, IB_mx_ptr, o);
}

void
merge_join_int32
(
    int32_T   *A,
    int32_T   *B,
    mxArray  **IA_mx_ptr,
    mxArray  **IB_mx_ptr,
    mwSize   A_size, 
    mwSize   B_size)
{
    mwSize i = 0;
    mwSize j = 0;
    mwSize o = 0;
    
    mwSize ii;
    mwSize jj;
    
    /* allocate output */
    mwSize o_sz = A_size + B_size;
    *IA_mx_ptr = mxCreateDoubleMatrix(o_sz, 1, mxREAL);
    *IB_mx_ptr = mxCreateDoubleMatrix(o_sz, 1, mxREAL);
    
    double *IA = mxGetPr(*IA_mx_ptr);
    double *IB = mxGetPr(*IB_mx_ptr);
    
    
    while(i != A_size && j != B_size) {
        if     (A[i] > B[j]) ++j;
        else if(A[i] < B[j]) ++i;
        else {
            o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, i, j, o, o_sz);
            ++o;
            jj = j+1;
            while(jj != B_size && A[i] == B[jj]) {
                o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, i, jj, o, o_sz);
                ++o;
                ++jj;
            }
            ii = i+1;
            while(ii != A_size && B[j] == A[ii]) {
                o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, ii, j, o, o_sz);
                ++o;
                ++ii;
            }
            ++i; ++j;
        }
    }
    
    /* resize the output to avoid trailing zeros */
    shrink_out(IA_mx_ptr, IB_mx_ptr, o);
}

void
merge_join_uint32
(
    uint32_T   *A,
    uint32_T   *B,
    mxArray  **IA_mx_ptr,
    mxArray  **IB_mx_ptr,
    mwSize   A_size, 
    mwSize   B_size)
{
    mwSize i = 0;
    mwSize j = 0;
    mwSize o = 0;
    
    mwSize ii;
    mwSize jj;
    
    /* allocate output */
    mwSize o_sz = A_size + B_size;
    *IA_mx_ptr = mxCreateDoubleMatrix(o_sz, 1, mxREAL);
    *IB_mx_ptr = mxCreateDoubleMatrix(o_sz, 1, mxREAL);
    
    double *IA = mxGetPr(*IA_mx_ptr);
    double *IB = mxGetPr(*IB_mx_ptr);
    
    
    while(i != A_size && j != B_size) {
        if     (A[i] > B[j]) ++j;
        else if(A[i] < B[j]) ++i;
        else {
            o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, i, j, o, o_sz);
            ++o;
            jj = j+1;
            while(jj != B_size && A[i] == B[jj]) {
                o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, i, jj, o, o_sz);
                ++o;
                ++jj;
            }
            ii = i+1;
            while(ii != A_size && B[j] == A[ii]) {
                o_sz = out_set(&IA, &IB, IA_mx_ptr, IB_mx_ptr, ii, j, o, o_sz);
                ++o;
                ++ii;
            }
            ++i; ++j;
        }
    }
    
    /* resize the output to avoid trailing zeros */
    shrink_out(IA_mx_ptr, IB_mx_ptr, o);
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if(nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:merge_join:nrhs","Two inputs required.");
    }
    if(nlhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:merge_join:nlhs","Two outputs required.");
    }
    size_t A_size;
    size_t B_size;
    double *IA;
    double *IB;
    /* get dimensions of the input matrix */
    A_size = mxGetM(prhs[0]);
    B_size = mxGetM(prhs[1]);
    /* get a pointer to the real data in the output matrix */
    mxArray **IA_mx_ptr = &plhs[0];
    mxArray **IB_mx_ptr = &plhs[1];
   /* check type of input*/
    if(mxIsUint32(prhs[0]) && mxIsUint32(prhs[1])) {
        uint32_T *A;               /* Nx1 input matrix */
        uint32_T *B;               /* Nx1 input matrix */
        /* create a pointer to the real data in the input matrix  */
        A = (uint32_T*)mxGetData(prhs[0]);
        B = (uint32_T*)mxGetData(prhs[1]);
        /* call the computational routine */
        merge_join_uint32(A, B, IA_mx_ptr, IB_mx_ptr, (mwSize)A_size, (mwSize)B_size);
    } else if(mxIsInt32(prhs[0]) && mxIsInt32(prhs[1])) {
        int32_T *A;               /* Nx1 input matrix */
        int32_T *B;               /* Nx1 input matrix */
        /* create a pointer to the real data in the input matrix  */
        A = (int32_T*)mxGetData(prhs[0]);
        B = (int32_T*)mxGetData(prhs[1]);
        /* call the computational routine */
        merge_join_int32(A, B, IA_mx_ptr, IB_mx_ptr, (mwSize)A_size, (mwSize)B_size);
    } else if(mxIsDouble(prhs[0]) && mxIsDouble(prhs[1])) {
        double *A;               /* Nx1 input matrix */
        double *B;               /* Nx1 input matrix */
        /* create a pointer to the real data in the input matrix  */
        A = mxGetPr(prhs[0]);
        B = mxGetPr(prhs[1]);
        /* call the computational routine */
        merge_join_double(A, B, IA_mx_ptr, IB_mx_ptr, (mwSize)A_size, (mwSize)B_size);
    } else {
        mexErrMsgIdAndTxt("MyToolbox:merge_join:nlhs","Only uint32, int32 and double inputs are supported.");
    }
}


// mwSize
// expand_out
// (
//     mxArray  **IA_mx_ptr,
//     mxArray  **IB_mx_ptr,
//     mwSize   o_sz
// )
// {
//     //mwSize new_sz = o_sz * 2;
//     mwSize new_sz = o_sz * 20;
//     double *ptrA    = mxGetPr(*IA_mx_ptr);
//     double *ptrB    = mxGetPr(*IB_mx_ptr);
//     void   *newptrA = mxRealloc(ptrA, sizeof(double) * new_sz);
//     void   *newptrB = mxRealloc(ptrB, sizeof(double) * new_sz);
//     mxSetPr(*IA_mx_ptr, newptrA);
//     mxSetPr(*IB_mx_ptr, newptrB);
//     mxSetM(*IA_mx_ptr, new_sz);
//     mxSetM(*IB_mx_ptr, new_sz);
// 
//     return new_sz;
// }
