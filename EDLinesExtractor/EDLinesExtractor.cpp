
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
	//Mat imageRight = getImageFromMxArray(prhs[1], CV_8UC1);

	LineDescriptor lineDesc;
	ScaleLines linesLeft;
	//ScaleLines linesRight;
	lineDesc.GetLineDescriptor(imageLeft, linesLeft);
	//lineDesc.GetLineDescriptor(imageRight, linesRight);

	/*PairwiseLineMatching lineMatch;
	std::vector<unsigned int> matchResult;
	lineMatch.LineMatching(linesLeft, linesRight, matchResult);*/

	Mat linesOutput(linesLeft.size(), 5, CV_64FC1);

	for (int i = 0; i < linesLeft.size(); i++) {
		linesOutput.at<double>(i, 0) = (double) i;
		linesOutput.at<double>(i, 1) = (double) linesLeft[i][0].startPointX;
		linesOutput.at<double>(i, 2) = (double) linesLeft[i][0].startPointY;
		linesOutput.at<double>(i, 3) = (double) linesLeft[i][0].endPointX;
		linesOutput.at<double>(i, 4) = (double) linesLeft[i][0].endPointY;
	}

	plhs[0] = mxCreateDoubleMatrix(linesOutput.rows, linesOutput.cols, mxREAL);
	double *data_out = (double *) mxGetData(plhs[0]);
	openCVMat2MatlabMat(linesOutput, data_out);
}