#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // -- MATLAB --
    // function [vn,VN_v] = normvec(v,jacMethod)
    //
    // % NORMVEC Normalize vector.
    // %   NORMVEC(V) is the unit length vector in the same direction and sense as
    // %   V. It is equal to V/norm(V).
    // %
    // %   [NV,NV_v] = NORMVEC(V) returns the Jacobian of the normed vector wrt V.
    // %
    // %   [NV,NV_v] = NORMVEC(V,method) , with method ~= 0, uses a scalar diagonal
    // %   matrix as Jacobian, meaning that the vector has been just scaled to
    // %   length one by a scalar factor.
    
    // -- INPUT ARGUMENTS --
    // -- v --
    double* data = (double*) mxGetData(prhs[0]);
    size_t nRows = mxGetM(prhs[0]);
    size_t nCols = mxGetN(prhs[0]);
    Mat v = matlabMat2OpenCVMat(data, nRows, nCols);
    // -- jacMethod --
    int jacMethod = 0;
    if (nrhs > 1) {
        double* p = (double*) mxGetData(prhs[1]);
        jacMethod = (int) p[0];
    }
    
    // -- FUNCTION START --
    InputArray vArray(v);
    double n2 = v.dot(vArray);
    double n  = sqrt(n2);
    
    Mat vn = v / n;
    Mat VN_v;
    
    if (nlhs > 1) {
        int s = v.rows;
        
        if (nrhs > 1 && jacMethod != 0) { // % use scalar method (approx)
            VN_v = Mat::eye(nRows, nRows, CV_64FC1) / n;
        } else {
            double n3 = n * n2;
            VN_v = (n2 * Mat::eye(nRows, nRows, CV_64FC1) - v * v.t()) / n3;
        }
    }
    
    // -- OUTPUT ARGUMENTS --
    // -- vn --
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
        double* vn_out = mxGetPr(plhs[0]);
        openCVMat2MatlabMat(vn, vn_out);
    }    
    // -- VN_v --
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(nRows, nRows, mxREAL);
        double* VN_v_out = mxGetPr(plhs[1]);
        openCVMat2MatlabMat(VN_v, VN_v_out);
    }
}