#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"
#include "../CInterface/MatrixRoi.hpp"
#include "../DataManagement/split_c.hpp"
#include "./Rotations/q2R_c.hpp"
#include "./Rotations/q2Pi_c.hpp"
#include "./Rotations/pi2pc_c.hpp"
#include "./frame_c.hpp"


using namespace cv;

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // -- MATLAB --
    // function F = updateFrame(F, norm)
    // 
    // % UPDATEFRAME Update frame structure.
    // %   F = UPDATEFRAME(F)  updates the matrices R, Pi and Pc and the
    // %   vectors t, q, it and iq in the frame structure F. The only required
    // %   field in F is F.x, the 7-vector containing translation and orientation
    // %   quaternion of the frame F:
    // %
    // %       F.x = [t;q] = [x y z a b c d]'
    // %
    // %   where q = [a b c d]' must already be a unit vector.
    // %
    // %   F = UPDATEFRAME(F, true) forces quaternion normalization.
    // %
    // %   As a result, the following fields are created/updated as follows:
    // %
    // %       F.t  = F.x(1:3);
    // %       F.q  = F.x(4:7);
    // %       F.R  = q2R(F.q);
    // %       F.Rt = F.R';
    // %       F.Pi = q2Pi(F.q);
    // %       F.Pc = pi2pc(F.Pi);
    // %       F.it = -F.R'*F.t;
    // %       F.iq = [1;-1;-1;-1].*F.q;
    // %
    // %   See also FRAME, Q2R, Q2PI, PI2PC, SPLITFRAME.
    
    // -- INPUT ARGUMENTS --
    // -- F (structure) --
    Frame* F                = new Frame;
    int nFieldsPerStruct    = mxGetNumberOfFields(prhs[0]);

	// get all field names
	const char **fieldNames;
	fieldNames = (const char**) mxCalloc(nFieldsPerStruct, sizeof(*fieldNames));
	for (int field = 0; field < nFieldsPerStruct; field++){
        fieldNames[field] = mxGetFieldNameByNumber(prhs[0], field);
    }
    
    const mxArray* data;
    for (int iField = 0; iField < nFieldsPerStruct; iField++) {
        
		// only search for field 'x' - all other field are calculated based on 'x'
        if (strcmp(fieldNames[iField], "x") == 0) {
			data = mxGetFieldByNumber(prhs[0], 0, iField);
			F->x = getOpenCVMatFromMxArray(data);
			break;
        }                 
    }

    // -- norm (boolean) --
    bool norm = false;
    if (nrhs > 1) {
        bool* p = (bool*) mxGetData(prhs[1]);
        norm = (bool) p[0];
    }

	// -- FUNCTION START --
	if (nrhs > 1 && norm == true) {
		printf("(!) updateFrame: You called this function with norm == true. This is currently not implemented (!)\n");
		// TODO: F.x(4:7)  = normvec(F.x(4:7)); --> normvec currently only available as mex-function
	}

	F->t	= getRoi(F->x, 0, 2, 0, 0);	
	F->q	= getRoi(F->x, 3, 6, 0, 0);
	q2R(F->q, F->R);
	F->Rt	= F->R.t();
	F->Pi	= q2Pi(F->q);
	F->Pc	= pi2pc(F->Pi);

	// -- OUTPUT ARGUMENTS --
	// create output structure 
	int dim[1] = {1};
	plhs[0] = mxCreateStructArray(1, dim, nFieldsPerStruct, fieldNames);
	mxFree((void *)fieldNames);

	// fill output structure with values from F
	mxArray *data_mx;
	double *data_out;

	for (int i = 0; i < nFieldsPerStruct; i++) {

		const char* name = mxGetFieldNameByNumber(plhs[0], i);

		if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "x") == 0) {
			data_mx = mxCreateDoubleMatrix(F->x.rows, F->x.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->x, data_out);
		} else if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "t") == 0) {
			data_mx = mxCreateDoubleMatrix(F->t.rows, F->t.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->t, data_out);
		} else if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "q") == 0) {
			data_mx = mxCreateDoubleMatrix(F->q.rows, F->q.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->q, data_out);
		} else if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "R") == 0) {
			data_mx = mxCreateDoubleMatrix(F->R.rows, F->R.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->R, data_out);
		} else if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "Rt") == 0) {
			data_mx = mxCreateDoubleMatrix(F->Rt.rows, F->Rt.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->Rt, data_out);
		} else if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "Pi") == 0) {
			data_mx = mxCreateDoubleMatrix(F->Pi.rows, F->Pi.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->Pi, data_out);
		} else if (strcmp(mxGetFieldNameByNumber(plhs[0], i), "Pc") == 0) {
			data_mx = mxCreateDoubleMatrix(F->Pc.rows, F->Pc.cols, mxREAL);
			data_out = (double*) mxGetData(data_mx);
			openCVMat2MatlabMat(F->Pc, data_out);
		} else // found field name that does not get updated in this function
			continue;

		mxSetFieldByNumber(plhs[0], 0, i, data_mx);
	}
}