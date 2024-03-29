function Rob = map2rob(Rob)
% MAP2ROB  Update Rob structure from the Map information.
%   ROB = MAP2ROB(ROB) updates the structure ROB to reflect the information
%   contained in the golbal map Map.
%
%   See also UPDATEFRAME.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

global Map

% normalize quaternion - mean and cov
% *************************************************************************
% C/C++ PORTIERUNG
% *************************************************************************
% -- MATLAB --
%[Map.x(Rob.frame.r(4:7)), NQ_q] = normvec(Map.x(Rob.frame.r(4:7)),1);
% -- C/C++ --
[Map.x(Rob.frame.r(4:7)), NQ_q] = normvec_c(Map.x(Rob.frame.r(4:7)),1);

% -- MATLAB --
%Map.P(Rob.frame.r(4:7),Map.used) = NQ_q*Map.P(Rob.frame.r(4:7),Map.used);
%Map.P(Map.used,Rob.frame.r(4:7)) = Map.P(Map.used,Rob.frame.r(4:7))*NQ_q';
% -- C/C++ --
Map.P(Rob.frame.r(4:7),Map.used) = MatrixMul(NQ_q, Map.P(Rob.frame.r(4:7), Map.used));
Map.P(Map.used,Rob.frame.r(4:7)) = MatrixMul(Map.P(Map.used,Rob.frame.r(4:7)), MatrixTranspose(NQ_q));

% means
Rob.state.x = Map.x(Rob.state.r);
Rob.frame.x = Map.x(Rob.frame.r);

% -- MATLAB --
%Rob.frame   = updateFrame(Rob.frame);
% -- C/C++ --
frame_temp   = updateFrame_c(Rob.frame);
% only assign the values that have been affected by updateFrame - the
% others might be inconsistent
Rob.frame.x  = frame_temp.x;
Rob.frame.t  = frame_temp.t;
Rob.frame.q  = frame_temp.q;
Rob.frame.R  = frame_temp.R;
Rob.frame.Rt = frame_temp.Rt;
Rob.frame.Pi = frame_temp.Pi;
Rob.frame.Pc = frame_temp.Pc;
% *************************************************************************
% C/C++ PORTIERUNG
% *************************************************************************


% covariances
% Rob.state.P = Map.P(Rob.state.r,Rob.state.r);
% Rob.frame.P = Map.P(Rob.frame.r,Rob.frame.r);



% ========== End of function - Start GPL license ==========


%   # START GPL LICENSE

%---------------------------------------------------------------------
%
%   This file is part of SLAMTB, a SLAM toolbox for Matlab.
%
%   SLAMTB is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   SLAMTB is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with SLAMTB.  If not, see <http://www.gnu.org/licenses/>.
%
%---------------------------------------------------------------------

%   SLAMTB is Copyright:
%   Copyright (c) 2008-2010, Joan Sola @ LAAS-CNRS,
%   Copyright (c) 2010-2013, Joan Sola,
%   Copyright (c) 2014-    , Joan Sola @ IRI-UPC-CSIC,
%   SLAMTB is Copyright 2009 
%   by Joan Sola, Teresa Vidal-Calleja, David Marquez and Jean Marie Codol
%   @ LAAS-CNRS.
%   See on top of this file for its particular copyright.

%   # END GPL LICENSE

