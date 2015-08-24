
#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "LineDescriptor.hpp"
#include "PairwiseLineMatching.hpp"

using namespace cv;
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	// ----- read input from matlab -------------------------------------------
	// single line to be matched
	double *p_currentLmk = (double*) mxGetData(prhs[0]);
	// all line observations from next image
	double *p_rawObs     = (double*) mxGetData(prhs[1]);
	// ------------------------------------------------------------------------

	// ----- restore current landmark object ----------------------------------
	ScaleLines currentLmk;
	currentLmk.clear();
	currentLmk.resize(2);
	OctaveSingleLine octave;

	octave.startPointX = p_currentLmk[0];
	octave.startPointY = p_currentLmk[1];
	octave.endPointX   = p_currentLmk[2];
	octave.endPointY   = p_currentLmk[3];
	octave.direction   = p_currentLmk[4];
	octave.salience	   = p_currentLmk[5];
	octave.lineLength  = p_currentLmk[6];
	octave.numOfPixels = p_currentLmk[7];

	octave.sPointInOctaveX = octave.startPointX;
	octave.sPointInOctaveY = octave.startPointY;
	octave.ePointInOctaveX = octave.endPointX;
	octave.ePointInOctaveY = octave.endPointY;
	octave.octaveCount = 0;
		
	currentLmk[0].push_back(octave);
	currentLmk[1].push_back(octave);

	// ----- restore previous observation objects -----------------------------
	ScaleLines rawObs;
	int nRows = mxGetM(prhs[1]);
	int nCols = mxGetN(prhs[1]);
	rawObs.clear();
	rawObs.resize(nRows);

	for (int row = 0; row < nRows; row++) {
		octave.startPointX = p_rawObs[row + 0 * nRows];
		octave.startPointY = p_rawObs[row + 1 * nRows];
		octave.endPointX   = p_rawObs[row + 2 * nRows];
		octave.endPointY   = p_rawObs[row + 3 * nRows];
		octave.direction   = p_rawObs[row + 4 * nRows];
		octave.salience    = p_rawObs[row + 5 * nRows];
		octave.lineLength  = p_rawObs[row + 6 * nRows];
		octave.numOfPixels = p_rawObs[row + 7 * nRows];

		octave.sPointInOctaveX = octave.startPointX;
		octave.sPointInOctaveY = octave.startPointY;
		octave.ePointInOctaveX = octave.endPointX;
		octave.ePointInOctaveY = octave.endPointY;
		octave.octaveCount = 0;

		rawObs[row].push_back(octave);
	}
	// ------------------------------------------------------------------------

	// ----- perform line matching --------------------------------------------
	PairwiseLineMatching lineMatcher;
	std::vector<pair<unsigned int, unsigned int>> matchResult;
	lineMatcher.LineMatching(currentLmk, rawObs, matchResult);
	// ------------------------------------------------------------------------

	// ----- line matching output to matlab -----------------------------------
	Mat matchingResultsOutput(matchResult.size(), 3, CV_64FC1);
	for (int i = 0; i < matchResult.size(); i++) {
		matchingResultsOutput.at<double>(i, 0) = (double) i;
		matchingResultsOutput.at<double>(i, 1) = (double) matchResult.at(i).first;
		matchingResultsOutput.at<double>(i, 2) = (double) matchResult.at(i).second;
	}

	plhs[0] = mxCreateDoubleMatrix(matchingResultsOutput.rows, matchingResultsOutput.cols, mxREAL);
	double *data_out = (double*) mxGetData(plhs[0]);
	openCVMat2MatlabMat(matchingResultsOutput, data_out);
	// ------------------------------------------------------------------------
}