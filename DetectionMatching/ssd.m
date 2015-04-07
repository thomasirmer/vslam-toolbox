function sc = ssd(I,J,SII,SJJ)

% SSD  Sum of Squared Differences coefficient
%   SSD(I,J) computes the SSD score of matrices I and J:
%
%     SSD = 1/N*sum((I-J).^2)
%     
%   where N = prod(size(I)) 
%   and   size(I) = size(J)
%
%   SSD(I,J,SII,SJJ) accepts useful intermediate results
%   that permit to speed up the calculations. These are:
%
%     SII = sum(sum(I.*I))
%     SJJ = sum(sum(J.*J))
%
%   See also PATCHCORR, ZNCC, CENSUS

% (c) 2005 Joan Sola

if size(I) ~= size(J)
    error ('Matrices must be the same size.')
else
    switch nargin
        case {1,2}
            SII = sum(sum(I.*I));
            SJJ = sum(sum(J.*J));
        case 3
            SJJ = sum(sum(J.*J));
    end
    
    SIJ = sum(sum(I.*J));
    
    N   = numel(I);
    
    sc  = sqrt((SII+SJJ-2*SIJ)/(N+eps));

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

