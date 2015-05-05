
#ifndef SPLIT_HPP
#define SPLIT_HPP

	void split_c(Mat A, double** &R) {

		/* -- function parameters --
		@param double ** & R:
			Reference to a 2D array
			(must be a reference because array is created within this function
			otherwise new double... wouldn't be accessible from outside this function)
		*/

		// -- MATLAB --
		// % SPLIT  Split vectors into scalars, or matrices into row vectors.
		// %   [s1,s2,...,sn] = SPLIT(V), with V a vector, returns all its components
		// %   in scalars s1 ... sn. It is an error if numel(V) < nargout.
		// %
		// %   [v1,...,vn] = SPLIT(M), with M a matrix, returns its rows as separate
		// %   vectors v1 ... vn. It is an error if size(M,2) < nargout.
		//
		// %   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

		// -- FUNCTION START --
		double **R_; // temporary

		if (A.rows == 1) {			// A is a column vector --> split into scalars
			R_ = new double*[A.cols]; // create array of pointers to double arrays
			for (int col = 0; col < A.cols; col++) {
				R_[col] = new double[1]; // each double array contains just the one column value
				R_[col][0] = A.at<double>(0, col);
			}
		} else if (A.cols == 1) {	// A is a row vector --> split into scalars
			R_ = new double*[A.rows]; // create array of pointers to double arrays
			for (int row = 0; row < A.rows; row++) {
				R_[row] = new double[1]; // each double array contains just the one row value
				R_[row][0] = A.at<double>(row, 0);
			}
		} else if (A.rows > 1 && A.cols > 1) { // A is a matrix --> split rows into vectors

			printf("split_c - 'A' is a matrix - not tested yet! Please verify results!!");

			R_ = new double*[A.rows]; // create array of pointers to double arrays
			for (int row = 0; row < A.rows; row++) {
				R_[row] = new double[A.cols]; // each double array contains the column values for the current row 
				for (int col = 0; col < A.cols; col++) {
					R_[row][col] = A.at<double>(row, col);
				}
			}
		}

		R = R_; // assign temporary data to reference
	}

#endif // SPLIT_HPP