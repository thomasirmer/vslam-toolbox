#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // -- INPUT ARGUMENTS --
    // -- Matrix A data --
    double* mat_a = (double*) mxGetData(prhs[0]);
    // -- Matrix A number of rows --
    size_t nRows_a = mxGetM(prhs[0]);
    // -- Matrix A number of columns --
    size_t nCols_a = mxGetN(prhs[0]);
    
    // -- Matrix B data --
    double* mat_b = (double*) mxGetData(prhs[1]);
    // -- Matrix B number of rows --
    size_t nRows_b = mxGetM(prhs[1]);
    // -- Matrix B number of columns --
    size_t nCols_b = mxGetN(prhs[1]);
    
    // -- FUNCTION START --
    // Check matrix dimensions
    // (A * B = C) with A: l-m matrix
    //                  B: m-n matrix
    //                  C: l-n matrix
    if (nCols_a != nRows_b) {
        printf("Error: Invalid matrix dimensions.\n");
        return;
    }
    
    // Create OpenCV data types
    Mat A = matlabMat2OpenCVMat(mat_a, nRows_a, nCols_a);
    Mat B = matlabMat2OpenCVMat(mat_b, nRows_b, nCols_b);
    Mat C = A * B;
    
    // -- OUTPUT ARGUMENTS --
    // -- Matrix C data --
    plhs[0] = mxCreateDoubleMatrix(nRows_a, nCols_b, mxREAL);
    double* mat_c = mxGetPr(plhs[0]);
    openCVMat2MatlabMat(C, mat_c);
}