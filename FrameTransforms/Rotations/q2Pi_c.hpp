#ifndef Q2PI_HPP
#define Q2PI_HPP

Mat q2Pi(Mat q) {
	// -- MATLAB --
	// % Q2PI  Pi matrix construction from quaternion.
	// %   PI = Q2PI(Q) Jacobian submatrix PI from quaternion
	// %
	// %   Given:  Q     = [a b c d]'  the attitude quaternion
	// %           W     = [p q r]'    the angular rates vector
	// %           OMEGA = W2OMEGA(W)  a skew symetric matrix 
	// %
	// %   The output matrix:
	// %
	// %                |-b -c -d |
	// %           PI = | a -d  c |  
	// %                | d  a -b |
	// %                |-c  b  a |  
	// % 
	// %   is the Jacobian of OMEGA*Q with respect to W
	// %
	// %   See also W2OMEGA, Q2R
    //
	// %   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

	// -- FUNCTION START --
	Mat Pi = Mat(4, 3, CV_64FC1);

	if (q.rows == 1) { // q is a column vector
		// row 0
		Pi.at<double>(0, 0) = -q.at<double>(0,1);
		Pi.at<double>(0, 1) = -q.at<double>(0,2);
		Pi.at<double>(0, 2) = -q.at<double>(0,3);
		// row 1
		Pi.at<double>(1, 0) =  q.at<double>(0,0);
		Pi.at<double>(1, 1) = -q.at<double>(0,3);
		Pi.at<double>(1, 2) =  q.at<double>(0,2);
		// row 2
		Pi.at<double>(2, 0) =  q.at<double>(0,3);
		Pi.at<double>(2, 1) =  q.at<double>(0,0);
		Pi.at<double>(2, 2) = -q.at<double>(0,1);
		// row 3
		Pi.at<double>(3, 0) = -q.at<double>(0,2);
		Pi.at<double>(3, 1) =  q.at<double>(0,1);
		Pi.at<double>(3, 2) =  q.at<double>(0,0);
	} else if (q.cols == 1) { // q is a row vector
		// row 0
		Pi.at<double>(0, 0) = -q.at<double>(1,0);
		Pi.at<double>(0, 1) = -q.at<double>(2,0);
		Pi.at<double>(0, 2) = -q.at<double>(3,0);
		// row 1
		Pi.at<double>(1, 0) =  q.at<double>(0,0);
		Pi.at<double>(1, 1) = -q.at<double>(3,0);
		Pi.at<double>(1, 2) =  q.at<double>(2,0);
		// row 2
		Pi.at<double>(2, 0) =  q.at<double>(3,0);
		Pi.at<double>(2, 1) =  q.at<double>(0,0);
		Pi.at<double>(2, 2) = -q.at<double>(1,0);
		// row 3
		Pi.at<double>(3, 0) = -q.at<double>(2,0);
		Pi.at<double>(3, 1) =  q.at<double>(1,0);
		Pi.at<double>(3, 2) =  q.at<double>(0,0);
	}

	return Pi;
}

#endif // Q2PI_HPP