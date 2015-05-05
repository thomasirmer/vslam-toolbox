#ifndef PI2PC_HPP
#define PI2PC_HPP

Mat pi2pc(Mat Pi) {
	// -- MATLAB --
	// % PI2PC  Pi matrix to conjugated Pi matrix conversion.
	// %   PC = PI2PC(PI) converts the matrix
	// %       PI = QUAT2PI(Q)
	// %   into 
	// %       PC = QUAT2PI(QC)
	// %   where QC = Q2QC(Q), the conjugated quaternion.
	// %
	// % See also Q2QC, Q2PI, TOFRAME
	// 
	// %   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

	// -- FUNCTION START --
	Mat Pc = Mat(4, 3, CV_64FC1);

	// column 0
	Pc.at<double>(0,0) = -Pi.at<double>(0,0);
	Pc.at<double>(1,0) =  Pi.at<double>(1,0);
	Pc.at<double>(2,0) = -Pi.at<double>(2,0);
	Pc.at<double>(3,0) = -Pi.at<double>(3,0);
	// column 1
	Pc.at<double>(0,1) = -Pi.at<double>(0,1);
	Pc.at<double>(1,1) = -Pi.at<double>(1,1);
	Pc.at<double>(2,1) =  Pi.at<double>(2,1);
	Pc.at<double>(3,1) = -Pi.at<double>(3,1);
	// column 2
	Pc.at<double>(0,2) = -Pi.at<double>(0,2);
	Pc.at<double>(1,2) = -Pi.at<double>(1,2);
	Pc.at<double>(2,2) = -Pi.at<double>(2,2);
	Pc.at<double>(3,2) =  Pi.at<double>(3,2);

	return Pc;
}

#endif // PI2PC_HPP