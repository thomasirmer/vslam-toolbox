function [newId, meas, exp, inn] = detectLin(lmkIds, raw, Sen)

% SIMDETECTLIN  Detect 2D segment in simulated raw data.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

apps  = raw.segments.app;

[newIds,newIdsIdx] = setdiff(apps,lmkIds);

if ~isempty(newIds)
    newId    = newIds(1);
    newIdx   = newIdsIdx(1);
    
    % --------------------------------------------------------
    % ----- ACTIVE-SEARCH FOR REAL IMAGES IMPLEMENTATION -----
        
    % project landmark to sensor grid (--> 6.1.3 The active search - 1b)
    coords = raw.segments.coord(:,newIdx);
    meanX = (coords(1) + coords(3)) / 2;
    meanY = (coords(2) + coords(4)) / 2;
    
    lmk.occupiedCells = cell(Sen.imGrid.numCells(1), Sen.imGrid.numCells(2));
    
    breakLoop = false;
    for row = 1:Sen.imGrid.numCells(1)
        for col = 1:Sen.imGrid.numCells(2)
            if (meanX > Sen.imGrid.imPatchIndices{row, col}.x1 && ...
                    meanX < Sen.imGrid.imPatchIndices{row, col}.x2 && ...
                    meanY > Sen.imGrid.imPatchIndices{row, col}.y1 && ...
                    meanY < Sen.imGrid.imPatchIndices{row, col}.y2)
                lmk.occupiedCells{row, col} = 'occupied';
                breakLoop = true;
                break;
            end
        end
        if (breakLoop)
            break;
        end
    end
    
    % find randomly selected 'not occupied' cell (--> 6.1.3 The active search - 1c)
    cellNotFound = true;
    while (cellNotFound)
        randX = ceil(rand() * Sen.imGrid.numCells(2));
        randY = ceil(rand() * Sen.imGrid.numCells(1));
        
        if (strcmp(lmk.occupiedCells{randY, randX}, 'occupied') == false)
            cellNotFound = false;
        end
    end
        
    % extract image from cell
    imgCoords = Sen.imGrid.imPatchIndices{randX, randY};
    subimage = raw.img(imgCoords.x1:imgCoords.x2, imgCoords.y1:imgCoords.y2);
    
    % detect harris point (--> 6.1.3 The active search - 1d)
    [harrisPoint, harrisPointStrength] = harris_strongest(subimage);
    
    % set measurement results
    meas.y = harrisPoint;                    % TODO: meas.y wird von anderen Funktionen als 4x4 erwartet ?!?!?
    meas.R = harrisPointStrength^2 * eye(2); % ??? eigentlich soll hier pixnoise verwendet werden ???
    
    exp.e = meas.y;
    exp.E = meas.R;
    
    % ----- ACTIVE-SEARCH FOR REAL IMAGES IMPLEMENTATION -----
    % --------------------------------------------------------
    
    y = raw.segments.coord(:,newIdx);
    R = blkdiag(Sen.par.pixCov, Sen.par.pixCov);
    
    inn.z  = [0;0];
    inn.Z  = meas.R;
    
    
else
    
    newId  = [];
    meas.y = [];
    meas.R = [];
    exp.e  = [];
    exp.E  = [];
    inn.z  = [];
    inn.Z  = [];
    
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

