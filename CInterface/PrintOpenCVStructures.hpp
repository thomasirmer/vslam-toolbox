
#ifndef PRINT_OPENCV_STRUCTURES_HPP
#define PRINT_OPENCV_STRUCTURES_HPP

#include <opencvmex.hpp>

using namespace cv;

void printMat(Mat matrix) {
	int nRows = matrix.rows;
	int nCols = matrix.cols;

	for (int row = 0; row < nRows; row++) {
		for (int col = 0; col < nCols; col++) {
			printf("%f ", matrix.at<double>(row, col));
		}
		printf("\n");
	}
}

#endif // PRINT_OPENCV_STRUCTURES_HPP