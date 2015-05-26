#include <opencvmex.hpp>

#include "../../CInterface/DataTypeParsing.hpp"
#include "../../CInterface/PrintOpenCVStructures.hpp"
#include "../../CInterface/MatrixRoi.hpp"

using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	// -- MATLAB --
	/*
	% V2Q Rotaiton vector to quaternion conversion.
	%   [Q,Qv] = V2Q(V) returns the quaternion Q correscponding to the rotation
	%   encoded in rotation vector V, and the associated Jacobian Qv = dQ/dV.
	*/

	// -- INPUT ARGUMENTS --
	Mat v = getOpenCVMatFromMxArray(prhs[0]);
}