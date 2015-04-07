function [l, L_rf, L_sf, L_k, L_seg, L_n] = ...
    retroProjIdpLinFromPinHoleOnRob(Rf, Sf, k, seg, n)

% RETROPROJIDPLINFROMPINHOLEONROB retroprj Idp Line from pinhole on robot.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.


if nargout == 1
    
    ls = invPinHoleIdpLin(k, seg, n) ;
    lr = fromFrameIdpLin(Sf, ls);
    l  = fromFrameIdpLin(Rf, lr);
    
else % Jacobians requested
    
    [ls, LS_seg, LS_n, LS_k] = invPinHoleIdpLin(seg, n, k) ;
    [lr, LR_sf, LR_ls]       = fromFrameIdpLin(Sf, ls);
    [l, L_rf, L_lr]          = fromFrameIdpLin(Rf, lr);

    L_sf  = L_lr*LR_sf;
    L_ls  = L_lr*LR_ls;
    L_k   = L_ls*LS_k;
    L_seg = L_ls*LS_seg;
    L_n   = L_ls*LS_n;
    
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

