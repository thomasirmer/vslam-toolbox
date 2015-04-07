function [l,L_k,L_hm,L_beta] = invPinHoleAPlucker(Sk,hm,beta)

% INVPINHOLEAPLUCKER Retro-projects anchored plucker line
%   INVPINHOLEAPLUCKER(K,L,BETA) retro-projects the anchored Plucker line
%   from the homogeneous 2D line HM and a pin hole camera K=[u0;v0;au;av]
%   at the origin. BETA specifies the unobservable direction of the line.
%   BETA is a 2-vector expressed in the plane base given by
%   PLANEVEC2PLANEBASE.
%
%   [l,L_k,L_hm,L_beta] = ... returns Jacobians wrt K, HM and BETA.
%
%   See also INVPINHOLEPLUCKER.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

if nargout == 1
    % Plucker line 
    pl = invPinHolePlucker(Sk,hm,beta) ;

    % anchored Plucker line L
    l = anchorPlucker(pl,[0;0;0]);

else
    % Plk lin in Sensor Frame
    [pl, PL_k, PL_hm, PL_beta] = invPinHolePlucker(Sk,hm,beta) ;
    
    % L in sensor frame
    [l, L_pl] = anchorPlucker(pl,[0;0;0]);
    
    % chain rule
    L_k    = L_pl * PL_k;
    L_hm   = L_pl * PL_hm;
    L_beta = L_pl * PL_beta;

end



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

