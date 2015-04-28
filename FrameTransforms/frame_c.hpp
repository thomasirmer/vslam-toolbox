
#ifndef FRAME_C_HPP
#define FRAME_C_HPP

#include <opencvmex.hpp>

// -- MATLAB --
// % FRAME  Help on frames for the FrameTransforms/ toolbox.
// %
// %   Frames are Matlab structures used to store data belonging to 3D frames.
// %   We do this to avoid having to compute multiple times rotation matrices
// %   and other frame-related constructions. 
// %
// %   A frame is specified with a translation vector and an orientation
// %   quaternion (see QUATERNION), giving a state vector of 7 components.
// %   This is the essential frame information, stored in field frame.x. From
// %   this 7-vector, all other fields are created or updated using the
// %   UPDATEFRAME function.
// %
// %   The fields of a frame structure are:
// %
// %     .x      * the state 7-vector
// %     .t      * translation vector,     t  = x(1:3)
// %     .q      * orientation quaternion, q  = x(4:7)
// %     .R      * rotation matrix,        R  = q2R(q)
// %     .Rt     * transposed R,           Rt = R'
// %     .it     * inverse position,       it = -Rt*t
// %     .iq     * inverse or conjugate quaternion, iq = q2qc(q)
// %     .Pi     * PI matrix,              Pi = q2Pi(q)
// %     .Pc     * conjugate PI matrix,    Pc = q2Pi(iq)
// % 
// %   Fields .Pi and .Pc are used to help computing some Jacobians of
// %   functions involving frame transformations. See TOFRAME in case you are
// %   interested in finding out more.
// %
// %   See also UPDATEFRAME, SPLITFRAME, QUATERNION, TOFRAME, FROMFRAME,
// %   COMPOSEFRAMES. 

struct Frame {
    Mat x, t, q, R, Rt, it, iq, Pi, Pc;
};

#endif // FRAME_C_HPP