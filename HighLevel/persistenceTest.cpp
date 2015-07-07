
#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

#include "../EDLinesExtractor/LineDescriptor.hpp"
#include "../EDLinesExtractor/PairwiseLineMatching.hpp"

static ScaleLines *lines = NULL;

void onExit()
{
	if (lines != NULL)
	{
		mexPrintf("Free memory...");
		mxFree(lines);
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	if (lines == NULL) // first call to this function
	{
		mexPrintf("First call! Allocating memory...\n");
		mexPrintf("Size of ScaleLines: %d bytes\n", sizeof(ScaleLines));
		lines = (ScaleLines*) mxMalloc(sizeof(ScaleLines));
		mexMakeMemoryPersistent(lines);
		mexAtExit(onExit);
	}
	else 
	{
		mexPrintf("Another call!\n");
	}
}