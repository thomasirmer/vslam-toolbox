
#ifndef MATRIX_ROI_HPP
#define MATRIX_ROI_HPP

#include <opencvmex.hpp>

#include "./PrintOpenCVStructures.hpp"

using namespace cv;

Mat getRoi(Mat matrix, int rowStart, int rowEnd, int colStart, int colEnd) {
	Mat roi(matrix, Rect(colStart, rowStart, colEnd - colStart + 1, rowEnd - rowStart + 1));
	return roi;
}

#endif // MATRIX_ROI_HPP