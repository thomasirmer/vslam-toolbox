
#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "./LineMatching/LineDescriptor.hpp"

using namespace cv;
using namespace std;

bool compareLength(LinesVec line1, LinesVec line2) {
	return (line1[0].lineLength > line2[0].lineLength);
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	// ----- read input from matlab ---------------------------------------------
	Mat image  = getImageFromMxArray(prhs[0], CV_8UC1);
	int nLines = (int) *(mxGetPr(prhs[1]));
	// ---------------------------------------------------------------------------

	// ----- detect lines --------------------------------------------------------
	LineDescriptor lineDesc;
	ScaleLines lineFeatures;
	lineDesc.GetLineDescriptor(image, lineFeatures);
	// ---------------------------------------------------------------------------

	// ----- take only the longest lines for performance reasons -----------------
	// sort lines descending based on length
	std::sort(lineFeatures.begin(), lineFeatures.end(), compareLength);

	// put n longest lines into new data structure
	ScaleLines nLongestLineFeatures;
	nLines = min(nLines, (int)lineFeatures.size());
	nLongestLineFeatures.resize(nLines);
	for (int i = 0; i < nLines; i++) {
		nLongestLineFeatures[i].push_back(lineFeatures[i][0]);
	}
	// ---------------------------------------------------------------------------

	// ----- load lines into output structure ------------------------------------
	Mat lineFeatureOutput(nLongestLineFeatures.size(), 8 + 72, CV_64FC1);
	for (int i = 0; i < nLongestLineFeatures.size(); i++) {
		lineFeatureOutput.at<double>(i, 0) = (double) nLongestLineFeatures[i][0].startPointX;
		lineFeatureOutput.at<double>(i, 1) = (double) nLongestLineFeatures[i][0].startPointY;
		lineFeatureOutput.at<double>(i, 2) = (double) nLongestLineFeatures[i][0].endPointX;
		lineFeatureOutput.at<double>(i, 3) = (double) nLongestLineFeatures[i][0].endPointY;
		lineFeatureOutput.at<double>(i, 4) = (double) nLongestLineFeatures[i][0].direction;
		lineFeatureOutput.at<double>(i, 5) = (double) nLongestLineFeatures[i][0].salience;
		lineFeatureOutput.at<double>(i, 6) = (double) nLongestLineFeatures[i][0].lineLength;
		lineFeatureOutput.at<double>(i, 7) = (double) nLongestLineFeatures[i][0].numOfPixels;

		// ----- line descriptor -------------------------------------------------
		for (int j = 0; j < 72; j++) {
			lineFeatureOutput.at<double>(i, j+8) = (double)nLongestLineFeatures[i][0].descriptor.at(j);
		}
	}
	// ---------------------------------------------------------------------------

	// ----- output to matlab ----------------------------------------------------
	double *data_out;

	plhs[0] = mxCreateDoubleMatrix(lineFeatureOutput.rows, lineFeatureOutput.cols, mxREAL);
	data_out = (double *)mxGetData(plhs[0]);
	openCVMat2MatlabMat(lineFeatureOutput, data_out);
	// ---------------------------------------------------------------------------
}