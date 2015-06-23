
#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

#include "LineDescriptor.hpp"
#include "PairwiseLineMatching.hpp"

using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	Mat imageLeft = getImageFromMxArray(prhs[0], CV_8UC1);
	Mat imageRight = getImageFromMxArray(prhs[1], CV_8UC1);

	LineDescriptor lineDesc;
	ScaleLines linesLeft;
	ScaleLines linesRight;
	lineDesc.GetLineDescriptor(imageLeft, linesLeft);
	lineDesc.GetLineDescriptor(imageRight, linesRight);

	PairwiseLineMatching lineMatch;
	std::vector<DMatch> matchResult;
	lineMatch.LineMatching(linesLeft, linesRight, matchResult);

	Mat linesOutputLeft(linesLeft.size(), 5, CV_64FC1);
	Mat linesOutputRight(linesRight.size(), 5, CV_64FC1);

	for (int i = 0; i < linesLeft.size(); i++) {
		linesOutputLeft.at<double>(i, 0) = (double) i;
		linesOutputLeft.at<double>(i, 1) = (double) linesLeft[i][0].startPointX;
		linesOutputLeft.at<double>(i, 2) = (double) linesLeft[i][0].startPointY;
		linesOutputLeft.at<double>(i, 3) = (double) linesLeft[i][0].endPointX;
		linesOutputLeft.at<double>(i, 4) = (double) linesLeft[i][0].endPointY;
	}

	for (int i = 0; i < linesLeft.size(); i++) {
		linesOutputRight.at<double>(i, 0) = (double) i;
		linesOutputRight.at<double>(i, 1) = (double) linesRight[i][0].startPointX;
		linesOutputRight.at<double>(i, 2) = (double) linesRight[i][0].startPointY;
		linesOutputRight.at<double>(i, 3) = (double) linesRight[i][0].endPointX;
		linesOutputRight.at<double>(i, 4) = (double) linesRight[i][0].endPointY;
	}

	double *data_out;

	plhs[0] = mxCreateDoubleMatrix(linesOutputLeft.rows, linesOutputLeft.cols, mxREAL);
	data_out = (double *) mxGetData(plhs[0]);
	openCVMat2MatlabMat(linesOutputLeft, data_out);

	plhs[1] = mxCreateDoubleMatrix(linesOutputRight.rows, linesOutputRight.cols, mxREAL);
	data_out = (double *) mxGetData(plhs[1]);
	openCVMat2MatlabMat(linesOutputRight, data_out);
}