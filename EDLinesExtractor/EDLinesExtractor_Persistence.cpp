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

Mat *image1 = NULL;

// ----- MEX INTERFACE FUNCTION -----
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	mexPrintf("set callback\n");
	mexAtExit(onExit);

	mexPrintf("get size\n");
	mwSize nBytes = sizeof(Mat);
	mexPrintf("size is %d bytes\n", nBytes);

	mexPrintf("allocate memory\n");
	image1 = (Mat*)mxMalloc(nBytes);

	if (image1->empty()) // first call
	{
		mexPrintf("copy data\n");
		copyDataToMat(*image1, prhs[0]);
		mexPrintf("memory persistence\n");
		mexMakeMemoryPersistent(image1);
	}
	else
	{
		mexPrintf("another call\n");
	}
}

// ----- FUNCTION DEFINITIONS -----
bool myCompare(LinesVec line1, LinesVec line2) {
	return (line1[0].lineLength > line2[0].lineLength);
}

void onExit() 
{
	mexPrintf("onExit\n");
	if (image1 != NULL) {
		mexPrintf("delete memory\n");
		delete image1;
	}
}