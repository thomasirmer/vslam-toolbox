
#ifndef Q2R_HPP
#define Q2R_HPP

#include "../../DataManagement/split_c.hpp"

void q2R(Mat q, Mat &R) {
	// -- MATLAB --
	// % Q2R  Quaternion to rotation matrix conversion.
	// %   R = Q2R(Q) builds the rotation matrix corresponding to the unit
	// %   quaternion Q. The obtained matrix R is such that the product:
	// %
	// %         rg = R * rb 
	// %
	// %   converts the body referenced vector  rb 
	// %     into the global referenced vector  rg
	// %
	// %   [R,Rq] = (...) returns the Jacobian wrt q.
	// %
	// %   See also QUATERNION, R2Q, Q2E, Q2V.
	//
	// %   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

	// -- FUNCTION START --

	// ---- split q into a,b,c,d ----------------------------------------------
	double **R_temp;
	split_c(q, R_temp);

	double a = R_temp[0][0];
	double b = R_temp[1][0];
	double c = R_temp[2][0];
	double d = R_temp[3][0];

	for ( int i = 0; i < 3; i++ )
		delete[] R_temp[i];
	delete[] R_temp;
	// ------------------------------------------------------------------------

	// ---- calculate R -------------------------------------------------------
	double aa = a * a;
	double ab = 2 * a * b;
	double ac = 2 * a * c;
	double ad = 2 * a * d;
	double bb = b * b;
	double bc = 2 * b * c;
	double bd = 2 * b * d;
	double cc = c * c;
	double cd = 2 * c * d;
	double dd = d * d;

	// R is actually a 2D array but in this case we handle it as 1D
	double *R_ = new double[9];

	// row 0
	R_[0] = aa + bb - cc - dd;
	R_[1] = bc - ad;
	R_[2] = bd + ac;
	// row 1
	R_[3] = bc + ad;
	R_[4] = aa - bb + cc -dd;
	R_[5] = cd - ab;
	// row 2
	R_[6] = bd - ac;
	R_[7] = cd + ab;
	R_[8] = aa - bb - cc + dd;

	R = Mat(3, 3, CV_64FC1, R_);
	// ------------------------------------------------------------------------
}

void q2R(Mat q, Mat &R, Mat &Rq) {
	// ---- version for this function with 2 output arguments (using Rq) ------

	q2R(q, R);

	// ---- split q into a,b,c,d ----------------------------------------------
	double **R_temp;
	split_c(q, R_temp);

	double a = R_temp[0][0];
	double b = R_temp[1][0];
	double c = R_temp[2][0];
	double d = R_temp[3][0];

	for ( int i = 0; i < 3; i++ )
		delete[] R_temp[i];
	delete[] R_temp;
	// ------------------------------------------------------------------------

	double a2 = 2 * a;
	double b2 = 2 * b;
	double c2 = 2 * c;
	double d2 = 2 * d;

	double *Rq_ = new double[36];

	// row 0
	Rq_[0]  = a2;
	Rq_[1]  = b2;
	Rq_[2]  = -c2;
	Rq_[3]  = -d2;
	// row 1
	Rq_[4]  = d2;
	Rq_[5]  = c2;
	Rq_[6]  = b2;
	Rq_[7]  = a2;
	// row 2
	Rq_[8]  = -c2;
	Rq_[9]  = d2;
	Rq_[10] = -a2;
	Rq_[11] = b2;
	// row 3
	Rq_[12] = -d2;
	Rq_[13] = c2;
	Rq_[14] = b2;
	Rq_[15] = -a2;
	// row 4
	Rq_[16] = a2;
	Rq_[17] = -b2;
	Rq_[18] = c2;
	Rq_[19] = -d2;
	// row 5
	Rq_[20] = b2;
	Rq_[21] = a2;
	Rq_[22] = d2;
	Rq_[23] = c2;
	// row 6
	Rq_[24] = c2;
	Rq_[25] = d2;
	Rq_[26] = a2;
	Rq_[27] = b2;
	// row 7
	Rq_[28] = -b2;
	Rq_[29] = -a2;
	Rq_[30] = d2;
	Rq_[31] = c2;
	// row 8
	Rq_[32] = a2;
	Rq_[33] = -b2;
	Rq_[34] = -c2;
	Rq_[35] = d2;

	Rq = Mat(9, 4, CV_64FC1, Rq_);
	// ------------------------------------------------------------------------
}

#endif // Q2R_HPP