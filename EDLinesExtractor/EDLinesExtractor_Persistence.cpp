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

// ----- MEX INTERFACE FUNCTION -----
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{

}

// ----- FUNCTION DEFINITIONS -----

bool myCompare(LinesVec line1, LinesVec line2) {
	return (line1[0].lineLength > line2[0].lineLength);
}

void onExit() 
{

}