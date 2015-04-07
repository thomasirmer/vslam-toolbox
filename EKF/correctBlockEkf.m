function correctBlockEkf(r,H,Inn)

% CORRECTBLOCKEKF  Correct in block-defined EKF.
%   CORRECTBLIOCKEKF(r,H,INN) performs a correction step to global map Map
%   by using the observation Jacobian H, referring to range r in the map,
%   and innovation INN. 
%
%   INN is a structure containing:
%       .z      the innovation,         z  = y-h(x)
%       .Z      its covariances matrix, Z  = HPH' + R
%       .iZ     the inverse covariance, iZ = Z^-1.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

global Map

% Kalman gain
K = Map.P(Map.used,r) * H' * Inn.iZ;   % K = PH'Z^-1

% mean and cov. updates
Map.x(Map.used)          = Map.x(Map.used)          + K*Inn.z;
Map.P(Map.used,Map.used) = Map.P(Map.used,Map.used) - K*Inn.Z*K';

% Force symmetry
%
% NOTE: this line of code has been moved to correctKnownLmks so that it is
% performed once per SLAM iteration, and not once per landmark correction.
% This is so done for speed reasons.
%
% Map.P(Map.used,Map.used) = (Map.P(Map.used,Map.used) + Map.P(Map.used,Map.used)')/2;



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

