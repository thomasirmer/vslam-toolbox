
#ifndef DATA_TYPE_PRASING_HPP
#define DATA_TYPE_PRASING_HPP

#include <opencvmex.hpp>

#include "./PrintOpenCVStructures.hpp"

using namespace cv;

// MATLAB --> OPENCV

Mat matlabMat2OpenCVMat(double* data, size_t nRows, size_t nCols) {
    
    // create Mat object as double matrix with given rows and columns
    Mat matrix(nRows, nCols, CV_64FC1);
    
    // fill matrix data with given data.
    for (int row = 0; row < nRows; row++) {
		for (int col = 0; col < nCols; col++) {
			matrix.at<double>(row, col) = data[col * nRows + row];
		}
	}
    
    return matrix;
}

Mat getOpenCVMatFromMxArray(const mxArray* mxArrayData) {
    double* data = (double*) mxGetData(mxArrayData);
    size_t nRows = mxGetM(mxArrayData);
    size_t nCols = mxGetN(mxArrayData);
    return matlabMat2OpenCVMat(data, nRows, nCols);
}

// OPENCV --> MATLAB

void openCVMat2MatlabMat(Mat matrix, double* data) {
    int nRows = matrix.rows;
    int nCols = matrix.cols;
    
    for (int row = 0; row < nRows; row++) {
		for (int col = 0; col < nCols; col++) {
			data[col * nRows + row] = matrix.at<double>(row, col);
		}
	}
}

#endif // DATA_TYPE_PRASING_HPP