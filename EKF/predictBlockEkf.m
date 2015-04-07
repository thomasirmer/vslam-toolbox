function predictBlockEkf(r,F_r,U,F_u)

% PREDICTBLOCKEKF  Covariance predict in block-defined EKF.
%   PREDICTBLIOCKEKF(r,F_r,U,F_u) performs a covariance prediction step to
%   global map Map by using the prediction Jacobian F_r, referring to range
%   r in the map, and perturbation covariances matrix U and Jacobian F_u.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

global Map

m = Map.used; % caution: the range 'm' includes the range 'r'

Map.P(r,m) = F_r * Map.P(r,m);
Map.P(m,r) = Map.P(m,r) * F_r';
Map.P(r,r) = Map.P(r,r) + F_u * U * F_u';



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

