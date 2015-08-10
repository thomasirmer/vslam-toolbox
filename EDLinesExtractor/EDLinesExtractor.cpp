
#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

#include "LineDescriptor.hpp"
#include "PairwiseLineMatching.hpp"
#include "LineMatchingChristopher\LineFeature.h"

using namespace cv;
using namespace std;

bool compareLength(LinesVec line1, LinesVec line2) {
	return (line1[0].lineLength > line2[0].lineLength);
}

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	// ----- read images from matlab ---------------------------------------------
	Mat imageLeft = getImageFromMxArray(prhs[0], CV_8UC1);
	Mat imageRight = getImageFromMxArray(prhs[1], CV_8UC1);
	// ---------------------------------------------------------------------------

	// ----- detect lines --------------------------------------------------------
	LineDescriptor lineDesc;
	ScaleLines linesLeft;
	ScaleLines linesRight;
	lineDesc.GetLineDescriptor(imageLeft, linesLeft);
	lineDesc.GetLineDescriptor(imageRight, linesRight);
	// ---------------------------------------------------------------------------

	// ----- take only the longest lines for performance reasons -----------------
	// sort lines descending based on length
	std::sort(linesLeft.begin(), linesLeft.end(), compareLength);
	std::sort(linesRight.begin(), linesRight.end(), compareLength);

	// put n longest lines into new data structure
	int n = 20;
	ScaleLines linesLeftLongest;
	ScaleLines linesRightLongest;
	linesLeftLongest.resize(n);
	linesRightLongest.resize(n);
	for (int i = 0; i < n; i++) {
		linesLeftLongest[i].push_back(linesLeft[i][0]);
		linesRightLongest[i].push_back(linesRight[i][0]);
	}
	// ---------------------------------------------------------------------------

	// ----- match lines ---------------------------------------------------------
	PairwiseLineMatching lineMatch;
	std::vector<pair<unsigned int, unsigned int>> matchResult;
	lineMatch.LineMatching(linesLeftLongest, linesRightLongest, matchResult);
	// ---------------------------------------------------------------------------

	// ----- load lines into output structure ------------------------------------
	Mat linesOutputLeft(linesLeftLongest.size(), 6, CV_64FC1);
	for (int i = 0; i < linesLeftLongest.size(); i++) {
		linesOutputLeft.at<double>(i, 0) = (double) i;
		linesOutputLeft.at<double>(i, 1) = (double) linesLeftLongest[i][0].startPointX;
		linesOutputLeft.at<double>(i, 2) = (double) linesLeftLongest[i][0].startPointY;
		linesOutputLeft.at<double>(i, 3) = (double) linesLeftLongest[i][0].endPointX;
		linesOutputLeft.at<double>(i, 4) = (double) linesLeftLongest[i][0].endPointY;
		linesOutputLeft.at<double>(i, 5) = (double) linesLeftLongest[i][0].lineLength;
	}

	Mat linesOutputRight(linesRightLongest.size(), 6, CV_64FC1);
	for (int i = 0; i < linesRightLongest.size(); i++) {
		linesOutputRight.at<double>(i, 0) = (double) i;
		linesOutputRight.at<double>(i, 1) = (double) linesRightLongest[i][0].startPointX;
		linesOutputRight.at<double>(i, 2) = (double) linesRightLongest[i][0].startPointY;
		linesOutputRight.at<double>(i, 3) = (double) linesRightLongest[i][0].endPointX;
		linesOutputRight.at<double>(i, 4) = (double) linesRightLongest[i][0].endPointY;
		linesOutputRight.at<double>(i, 5) = (double) linesRightLongest[i][0].lineLength;
	}

	Mat matchingResultsOutput(matchResult.size(), 3, CV_64FC1);
	for (int i = 0; i < matchResult.size(); i++) {
		matchingResultsOutput.at<double>(i, 0) = (double)i;
		matchingResultsOutput.at<double>(i, 1) = (double)matchResult.at(i).first;
		matchingResultsOutput.at<double>(i, 2) = (double)matchResult.at(i).second;
	}

	/*Mat matchingResultsOutput(matchResult.size(), 5, CV_64FC1);
	for (int i = 0; i < matchResult.size(); i++) {
		matchingResultsOutput.at<double>(i, 0) = (double) i;
		matchingResultsOutput.at<double>(i, 1) = (double) matchResult.at(i).distance;
		matchingResultsOutput.at<double>(i, 2) = (double) matchResult.at(i).imgIdx;
		matchingResultsOutput.at<double>(i, 3) = (double) matchResult.at(i).trainIdx;
		matchingResultsOutput.at<double>(i, 4) = (double) matchResult.at(i).queryIdx;
	}*/
	// ---------------------------------------------------------------------------

	// ----- output to matlab ----------------------------------------------------
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
	// ---------------------------------------------------------------------------
}