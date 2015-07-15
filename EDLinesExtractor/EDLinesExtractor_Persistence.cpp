#include <opencvmex.hpp>

#include <math.h>	
#include <time.h>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

#include "LineDescriptor.hpp"
#include "PairwiseLineMatching.hpp"

using namespace cv;
using namespace std;

// ----- FUNCTION DECLARATIONS -----
bool myCompare(LinesVec, LinesVec);
void onExit();

// ----- MEMORY VARIABLES -----
static mxArray *addrImage1;
static mxArray *addrImage2;

// ----- MEX INTERFACE FUNCTION -----
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{

	mexPrintf("");
	
	return;
	
	
	/*

	Mat image1;
	Mat image2;

	if (image1.empty() && image2.empty()) // first call
	{
		mexPrintf("first call!");
		
		image1 = getImageFromMxArray(prhs[0], CV_8UC1);
		
		mexMakeMemoryPersistent(&image1);
		
	} 
	else if (image2.empty()) // second call
	{
		mexPrintf("second call!");
		
		image2 = getImageFromMxArray(prhs[0], CV_8UC1);
		mexMakeMemoryPersistent(&image1);
		
	}
	else // all other calls
	{
		mexPrintf("another call!");
		
		image1 = image2;
		image2 = getImageFromMxArray(prhs[0], CV_8UC1);
		
	}
	
	*/

}

// ----- FUNCTION DEFINITIONS -----
bool myCompare(LinesVec line1, LinesVec line2) {
	return (line1[0].lineLength > line2[0].lineLength);
}

void onExit() 
{

}