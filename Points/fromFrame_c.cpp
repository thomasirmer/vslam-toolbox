#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"
#include "../CInterface/MatrixRoi.hpp"

#include "../FrameTransforms/frame_c.hpp"

using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	// -- MATLAB --
	/*
	% FROMFRAME  Express in global frame a set of points from a local frame.
	%   FROMFRAME(F,P_F)  takes the F-referenced points matrix P_F and 
	%   returns it in the frame where F is referenced.
	%   P_F is a points matrix defined as 
	%     P_F = [P1 P2 ... PN], where
	%     Pi  = [xi;yi;zi]
	%
	%   F is either a structure containing at least:
	%     t : frame position
	%     q : frame orientation quaternion
	%     R : rotation matrix
	%     Pi: Pi matrix
	%
	%   or a 7-vector F = [t;q].
	%
	%   [p_W,Ff,Fp] = ... returns the Jacobians of fromFrame:
	%     Ff: wrt the frame [t;q]
	%     Fp: wrt the point P_F
	%   Note that this is only available for single points.
	%
	%   See also FRAME, TOFRAME, Q2PI, QUATERNION, UPDATEFRAME, SPLITFRAME.

	%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.


	%   [1] Joan Sola, "Towards visual localization, mapping and moving objects
	%   tracking by a moible robot," PhD dissertation, pages 181-183, Institut
	%   National Politechnique de Toulouse, 2007.
	*/

	// -- INPUT ARGUMENTS --
	// -- F (structure) --
	Frame* F			= new Frame;
	int nFieldPerStruct = mxGetNumberOfFields(prhs[0]);

	// get all field names
	const char **fieldNames;
	fieldNames = (const char **) mxCalloc(nFieldPerStruct, sizeof(*fieldNames));
	for (int field = 0; field < nFieldPerStruct; field++) {
		fieldNames[field] = mxGetFieldNameByNumber(prhs[0], field);
	}

	// get all fields
	const mxArray *data;
	for (int field = 0; field < nFieldPerStruct; field++) {
		data = mxGetFieldByNumber(prhs[0], 0, field);

		if (strcmp(fieldNames[field], "t") == 0) {
			F->t = getOpenCVMatFromMxArray(data);
		} else if (strcmp(fieldNames[field], "q") == 0) {
			F->q = getOpenCVMatFromMxArray(data);
		} else if (strcmp(fieldNames[field], "R") == 0) {
			F->R = getOpenCVMatFromMxArray(data);
		} else if (strcmp(fieldNames[field], "Rt") == 0) {
			F->Rt = getOpenCVMatFromMxArray(data);
		} else if (strcmp(fieldNames[field], "Pi") == 0) {
			F->Pi = getOpenCVMatFromMxArray(data);
		}
	}

	// -- p_F (vector) --
	Mat p_F = getOpenCVMatFromMxArray(prhs[1]);

	// -- FUNCTION START --
	int size = p_F.cols;

	Mat p_W;
	Mat Ff;
	Mat Fp;

	if (size == 1) { // % one point
		p_W = F->R * p_F + F->t;

		if (nlhs > 1) { // % Jacobians. See [1] for details.
			Mat s = 2 * F->Pi * p_F;

			Mat Ft = Mat::eye(3, 3, CV_64FC1);

			Mat Fq = Mat(3, 4, CV_64FC1);
			// row 0
			Fq.at<double>(0,0) =  s.at<double>(0,1);
			Fq.at<double>(0,1) = -s.at<double>(0,0);
			Fq.at<double>(0,2) =  s.at<double>(0,3);
			Fq.at<double>(0,3) = -s.at<double>(0,2);
			// row 1
			Fq.at<double>(1,0) =  s.at<double>(0,2);
			Fq.at<double>(1,1) = -s.at<double>(0,3);
			Fq.at<double>(1,2) = -s.at<double>(0,0);
			Fq.at<double>(1,3) =  s.at<double>(0,1);
			// row 2
			Fq.at<double>(2,0) =  s.at<double>(0,3);
			Fq.at<double>(2,1) =  s.at<double>(0,2);
			Fq.at<double>(2,2) = -s.at<double>(0,1);
			Fq.at<double>(2,3) = -s.at<double>(0,0);

			Fp = F->R;

			hconcat(Ft, Fq, Ff);
		} 
	} else { // % multiple points
		p_W = F->R * p_F + repeat(F->t, 1, size);
	}

	// -- OUTPUT ARGUMENTS --
	// -- p_W --
	plhs[0] = mxCreateDoubleMatrix(p_W.rows, p_W.cols, mxREAL);
	double *p_W_out = (double *) mxGetData(plhs[0]);
	openCVMat2MatlabMat(p_W, p_W_out);

	// -- Ff --
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(Ff.rows, Ff.cols, mxREAL);
		double *Ff_out = (double *) mxGetData(plhs[1]);
		openCVMat2MatlabMat(Ff, Ff_out);
	}

	// -- Fp --
	if (nlhs > 1) {
		plhs[2] = mxCreateDoubleMatrix(Fp.rows, Fp.cols, mxREAL);
		double *Fp_out = (double *) mxGetData(plhs[2]);
		openCVMat2MatlabMat(Fp, Fp_out);
	}
}