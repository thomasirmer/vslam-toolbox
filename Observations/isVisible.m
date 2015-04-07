function vis = isVisible(pix,depth,imSize,mrg)

% ISVISIBLE  Points visible from pinHole camera.
%   ISVISIBLE(PIX,DEPTH,IMSIZE) returns TRUE for those pixels in pixels
%   matrix PIX that are in front of the camera (DEPTH>0) and within the
%   limits defined by IMSIZE.
%  
%   ISVISIBLE(PIX,DEPTH,IMSIZE,MARGIN) considers a point inside the image
%   if it is further than MARGIN units from its borders.
%  
%   IMSIZE is defined by IMSIZE = [hsize vsize].
%  
%   PIX is a 2D-pixels matrix : PIX = [P1 ... PN] ,
%   with Pi = [ui;vi].
%
%   See also INSQUARE.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

%   TODO : check a better visibility criterion for negative inverse depth
%   management. Use the variance in RHO for this.


if nargin < 4
    mrg = 0;
end

vis = inSquare(pix,[0 imSize(1) 0 imSize(2)],mrg) & (depth > 0);



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

