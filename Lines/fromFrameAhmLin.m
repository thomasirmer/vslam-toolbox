function [l,L_f,L_lf] = fromFrameAhmLin(F,lf)

% FROMFRAMEAHMLIN  Transforms AHM line from local frame to global frame.
%   I = FROMFRAMEAHMLIN(F,LF) transforms the Anchored homogeneous line LF from the
%   local frame F to the global frame. The frame F must be specified via a
%   structure containing at least the fields F.t, F.q, F.R and F.Rt
%   (translation, quaternion, rotation matrix and its transpose).
%
%   [I,L_f,L_lf] = FROMFRAMEAHMLIN(...) returns the Jacobians wrt F and IF.

%   Copyright 2009 Teresa Vidal.

if nargout == 1
    [p1f,p2f] = ahmLin2ahmPnts(lf);

    p1 = fromFrameAhm(F,p1f);
    p2 = fromFrameAhm(F,p2f);

    l  = [p1;p2(4:7,:)];

else
    [p1f,p2f,P1f_lf,P2f_lf] = ahmLin2ahmPnts(lf);

    [p1,P1_f,P1_p1f] = fromFrameAhm(F,p1f);
    [p2,P2_f,P2_p2f] = fromFrameAhm(F,p2f);

    l  = [p1;p2(4:7)];
    
    % Jacobians
    L_f   = [P1_f;P2_f(4:7,:)];
    P2_lf = P2_p2f*P2f_lf;
    L_lf  = [P1_p1f*P1f_lf;P2_lf(4:7,:)];

    
end

return

%% jac

syms x y z a b c d X Y Z A1 B1 C1 R1 A2 B2 C2 R2 real
F.x = [x;y;z;a;b;c;d];
F   = updateFrame(F);
l_F = [X;Y;Z;A1;B1;C1;R1;A2;B2;C2;R2];

[l,L_f,L_lf] = fromFrameAhmLin(F,l_F);

simplify(L_f  - jacobian(l,F.x))
simplify(L_lf - jacobian(l,l_F))



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

