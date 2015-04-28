#include <opencvmex.hpp>

#include "../CInterface/DataTypeParsing.hpp"
#include "../CInterface/PrintOpenCVStructures.hpp"

#include "./frame_c.hpp";

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
    mwSize nStructElements  = mxGetNumberOfElements(prhs[0]);
    int nFieldsPerStruct    = mxGetNumberOfFields(prhs[0]);
    //mxClassID* classIDFlags = (mxClassID*)   mxCalloc(nFieldsPerStruct, sizeof(mxClassID));
    //const char** fieldNames = (const char**) mxCalloc(nFieldsPerStruct, sizeof(*fieldNames));
    
    const mxArray* data;
    for (int iField = 0; iField < nFieldsPerStruct; iField++) {
        
        data = mxGetFieldByNumber(prhs[0], 0, iField);
        const char* fieldName = mxGetFieldNameByNumber(prhs[0], iField);
        
        if (strcmp(fieldName, "x") == 0) {
            printf("field %d: %s\n", iField, fieldName);
            Mat x_ = getOpenCVMatFromMxArray(data);
            F->x = x_;
            printMat(F->x);
        }
        
        delete fieldName;
        
        //classIDFlags[iField] = mxGetClassID(data);
    }
    
    // -- norm (boolean) --
    bool norm = false;
    if (nrhs > 1) {
        bool* p = (bool*) mxGetData(prhs[1]);
        norm = (bool) p[0];
    }
}