
#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "LineDescriptor.hpp"

using namespace cv;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	// ----- read input from matlab ---------------------------------------------
	// single line to be matched
	double *p_currentLmk = (double*) mxGetData(prhs[0]);
	// all line observations from next image
	//double *p_rawObs     = (double*) mxGetData(prhs[1]);
	// ---------------------------------------------------------------------------

	// ----- restore Scalelines object -------------------------------------------
	ScaleLines currentLmk;
	OctaveSingleLine octave;
	octave.startPointX = p_currentLmk[0];
	octave.startPointY = p_currentLmk[1];
	octave.endPointX   = p_currentLmk[2];
	octave.endPointY   = p_currentLmk[3];
	// ---------------------------------------------------------------------------

	//// ----- take only the longest lines for performance reasons -----------------
	//// sort lines descending based on length
	//std::sort(lineFeatures.begin(), lineFeatures.end(), compareLength);

	//// put n longest lines into new data structure
	//ScaleLines nLongestLineFeatures;
	//nLongestLineFeatures.resize(nLines);
	//for (int i = 0; i < nLines; i++) {
	//	nLongestLineFeatures[i].push_back(lineFeatures[i][0]);
	//}
	//// ---------------------------------------------------------------------------

	//// ----- load lines into output structure ------------------------------------
	//Mat lineFeatureOutput(nLongestLineFeatures.size(), 5, CV_64FC1);
	//for (int i = 0; i < nLongestLineFeatures.size(); i++) {
	//	lineFeatureOutput.at<double>(i, 0) = (double) nLongestLineFeatures[i][0].startPointX;
	//	lineFeatureOutput.at<double>(i, 1) = (double) nLongestLineFeatures[i][0].startPointY;
	//	lineFeatureOutput.at<double>(i, 2) = (double) nLongestLineFeatures[i][0].endPointX;
	//	lineFeatureOutput.at<double>(i, 3) = (double) nLongestLineFeatures[i][0].endPointY;
	//	lineFeatureOutput.at<double>(i, 4) = (double) nLongestLineFeatures[i][0].lineLength;
	//}
	//// ---------------------------------------------------------------------------

	//// ----- output to matlab ----------------------------------------------------
	//double *data_out;

	//plhs[0] = mxCreateDoubleMatrix(lineFeatureOutput.rows, lineFeatureOutput.cols, mxREAL);
	//data_out = (double *)mxGetData(plhs[0]);
	//openCVMat2MatlabMat(lineFeatureOutput, data_out);
	//// ---------------------------------------------------------------------------
}