#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"
#include "../CInterface/MatrixRoi.hpp"

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
    
	// -- roi dimensions --
	double* p;
    // -- from row --
	p = (double*) mxGetData(prhs[1]);
	int fromRow = (int) p[0];
	// -- to row --
	p = (double*) mxGetData(prhs[2]);
	int toRow = (int) p[0];
	// -- from column --
	p = (double*) mxGetData(prhs[3]);
	int fromCol = (int) p[0];
	// -- to row --
	p = (double*) mxGetData(prhs[4]);
	int toCol = (int) p[0];
    
    // -- FUNCTION START --
    // Create OpenCV data types
    Mat A = matlabMat2OpenCVMat(mat_in, nRows, nCols);


	// Get ROI
	Mat roi = getRoi(A, fromRow - 1, toRow - 1, fromCol - 1, toCol - 1);

	// -- OUTPUT ARGUMENTS --
	plhs[0] = mxCreateDoubleMatrix(roi.rows, roi.cols, mxREAL);
	double* roi_mat = mxGetPr(plhs[0]);
	openCVMat2MatlabMat(roi, roi_mat);
}