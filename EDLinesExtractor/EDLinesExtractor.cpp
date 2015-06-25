
#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

#include "LineDescriptor.hpp"
#include "PairwiseLineMatching.hpp"

using namespace cv;

bool myCompare(LinesVec line1, LinesVec line2) {
	return (line1[0].lineLength > line2[0].lineLength);
}

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	// ----- INPUT -----
	Mat imageLeft = getImageFromMxArray(prhs[0], CV_8UC1);
	Mat imageRight = getImageFromMxArray(prhs[1], CV_8UC1);

	// ----- FUNCTION START -----
	LineDescriptor lineDesc;
	ScaleLines lines1;
	ScaleLines lines2;
	lineDesc.GetLineDescriptor(imageLeft, lines1);
	lineDesc.GetLineDescriptor(imageRight, lines2);

	// ----- DEBUG -----
	// take only the n longest lines for performance reasons
	int n = 15;

	// sort lines descending based on length
	std::sort(lines1.begin(), lines1.end(), myCompare);
	std::sort(lines2.begin(), lines2.end(), myCompare);

	ScaleLines linesLeft;
	ScaleLines linesRight;
	linesLeft.resize(n);
	linesRight.resize(n);

	for (int i = 0; i < n; i++) {
		linesLeft[i].push_back(lines1[i][0]);
		linesRight[i].push_back(lines2[i][0]);
	}
	// ----- DEBUG -----

	PairwiseLineMatching lineMatch;
	std::vector<DMatch> matchResult;
	lineMatch.LineMatching(linesLeft, linesRight, matchResult);

	Mat linesOutputLeft(linesLeft.size(), 6, CV_64FC1);
	for (int i = 0; i < linesLeft.size(); i++) {
		linesOutputLeft.at<double>(i, 0) = (double) i;
		linesOutputLeft.at<double>(i, 1) = (double) linesLeft[i][0].startPointX;
		linesOutputLeft.at<double>(i, 2) = (double) linesLeft[i][0].startPointY;
		linesOutputLeft.at<double>(i, 3) = (double) linesLeft[i][0].endPointX;
		linesOutputLeft.at<double>(i, 4) = (double) linesLeft[i][0].endPointY;
		linesOutputLeft.at<double>(i, 5) = (double) linesLeft[i][0].lineLength;
	}

	Mat linesOutputRight(linesRight.size(), 6, CV_64FC1);
	for (int i = 0; i < linesRight.size(); i++) {
		linesOutputRight.at<double>(i, 0) = (double) i;
		linesOutputRight.at<double>(i, 1) = (double) linesRight[i][0].startPointX;
		linesOutputRight.at<double>(i, 2) = (double) linesRight[i][0].startPointY;
		linesOutputRight.at<double>(i, 3) = (double) linesRight[i][0].endPointX;
		linesOutputRight.at<double>(i, 4) = (double) linesRight[i][0].endPointY;
		linesOutputRight.at<double>(i, 5) = (double) linesRight[i][0].lineLength;
	}

	Mat matchingResultsOutput(matchResult.size(), 5, CV_64FC1);
	for (int i = 0; i < matchResult.size(); i++) {
		matchingResultsOutput.at<double>(i, 0) = (double) i;
		matchingResultsOutput.at<double>(i, 1) = (double) matchResult.at(i).distance;
		matchingResultsOutput.at<double>(i, 2) = (double) matchResult.at(i).imgIdx;
		matchingResultsOutput.at<double>(i, 3) = (double) matchResult.at(i).trainIdx;
		matchingResultsOutput.at<double>(i, 4) = (double) matchResult.at(i).queryIdx;
	}

	// ----- OUTPUT -----
	double *data_out;

	plhs[0] = mxCreateDoubleMatrix(linesOutputLeft.rows, linesOutputLeft.cols, mxREAL);
	data_out = (double *) mxGetData(plhs[0]);
	openCVMat2MatlabMat(linesOutputLeft, data_out);

	plhs[1] = mxCreateDoubleMatrix(linesOutputRight.rows, linesOutputRight.cols, mxREAL);
	data_out = (double *) mxGetData(plhs[1]);
	openCVMat2MatlabMat(linesOutputRight, data_out);

	plhs[2] = mxCreateDoubleMatrix(matchingResultsOutput.rows, matchingResultsOutput.cols, mxREAL);
	data_out = (double *) mxGetData(plhs[2]);
	openCVMat2MatlabMat(matchingResultsOutput, data_out);
}