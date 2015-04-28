#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{  
    // -- INPUT ARGUMENTS --
    // -- Matrix data --
    double* mat_in = (double*) mxGetData(prhs[0]);
    // -- Matrix A number of rows --
    size_t nRows = mxGetM(prhs[0]);
    // -- Matrix A number of columns --
    size_t nCols = mxGetN(prhs[0]);
    
    // -- FUNCTION START --
    // Create OpenCV data types
    Mat A = matlabMat2OpenCVMat(mat_in, nRows, nCols);
    
    // Transpose matrix
    A = A.t();
    
    // -- OUTPUT ARGUMENTS --
    // -- Matrix data --
    plhs[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
    double* mat_out = mxGetPr(plhs[0]);
    openCVMat2MatlabMat(A, mat_out);
}