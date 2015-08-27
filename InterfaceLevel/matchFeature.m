function Obs = matchFeature(Sen,Raw,Obs,Rob)

% MATCHFEATURE  Match feature.
% 	Obs = MATCHFEATURE(Sen,Raw,Obs) matches one feature in Raw to the predicted
% 	feature in Obs.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.

switch Obs.ltype(4:6)
    case 'Pnt'
        rawDataLmks = Raw.data.points;
        R = Sen.par.pixCov;
    case 'Lin'
        rawDataLmks = Raw.data.segments;
        R = blkdiag(Sen.par.pixCov,Sen.par.pixCov);
    otherwise
        error('??? Unknown landmark type ''%s''.',Obs.ltype);
end

switch Raw.type
    
    case {'simu','dump'}
        
        id  = Obs.lid;
        idx = find(rawDataLmks.app==id);
        
        if ~isempty(idx)
            Obs.meas.y   = rawDataLmks.coord(:,idx);
            Obs.meas.R   = R;
            Obs.measured = true;
            Obs.matched  = true;
        else
            Obs.meas.y   = zeros(size(Obs.meas.y));
            Obs.meas.R   = R;
            Obs.measured = false;
            Obs.matched  = false;
        end
        
    case 'image'
        % --------------------------------------------------------
        % ----- ACTIVE-SEARCH FOR REAL IMAGES IMPLEMENTATION -----
        id  = Obs.lid;
        idx = find(rawDataLmks.app==id);
        
        if ~isempty(idx)
            
            % --------------------------------------------------------
            % homogenous coordinates observation
            obs_l = Obs.exp.e(1);
            obs_m = Obs.exp.e(2);
            obs_n = Obs.exp.e(3);
            
            % --------------------------------------------------------
            % homogenous coordinates raw lines
            nLines = size(rawDataLmks.lines, 1);
            for i = 1:nLines
                [rawDataLmks.hmgCoords(:,i), ~] = seg2hmgLin(rawDataLmks.coord(:,i));
            end
            
            % --------------------------------------------------------
            % check all possible lines
            j = 1;
            possibleMatches = zeros(nLines, size(rawDataLmks.lines, 2) + 1);
            for i = 1:nLines
                % normalize homogenous coordinates of raw line to observation
                raw_l = rawDataLmks.hmgCoords(1,i);
                raw_m = rawDataLmks.hmgCoords(2,i);
                raw_n = rawDataLmks.hmgCoords(3,i);
                                
                norm = obs_n / raw_n;
                
                raw_l = raw_l * norm;
                raw_m = raw_m * norm;
                raw_n = raw_n * norm;
                
                % calculate 3-sigma neighborhood
                l_3sigma = 10 * sqrt(Obs.exp.E(1,1));
                m_3sigma = 10 * sqrt(Obs.exp.E(2,2));
                n_3sigma = 10 * sqrt(Obs.exp.E(3,3));
                
                % check if raw line in 3-sigma range
                if (raw_l > obs_l - l_3sigma && raw_l < obs_l + l_3sigma && ...
                    raw_m > obs_m - m_3sigma && raw_m < obs_m + m_3sigma && ...
                    raw_n > obs_n - n_3sigma && raw_n < obs_n + n_3sigma)
                    
                    possibleMatches(j,1:8) = rawDataLmks.lines(i,:);
                    possibleMatches(j,9) = i;
                    %plotSegLine(rawDataLmks.coord(:,i));
                    j = j + 1;
                end
            end
            
            %plotHmgLine(Obs.meas, Obs.exp);
                        
            matchedLine = MatchLines(Obs.meas.line, rawDataLmks.lines);
            
            Obs.meas.y   = possibleMatches(matchedLine(1,3)+1, 1:4)';
            Obs.meas.R   = R;
            Obs.measured = true;
            Obs.matched  = true;
        else
            Obs.meas.y   = zeros(size(Obs.meas.y));
            Obs.meas.R   = R;
            Obs.measured = false;
            Obs.matched  = false;
        end
%         end
%         id  = Obs.lid;
%         idx = find(rawDataLmks.app==id);
%         
%         if ~isempty(idx)
%             Obs.meas.y   = rawDataLmks.coord(:,idx);
%             Obs.meas.R   = R;
%             Obs.measured = true;
%             Obs.matched  = true;
%         else
%             Obs.meas.y   = zeros(size(Obs.meas.y));
%             Obs.meas.R   = R;
%             Obs.measured = false;
%             Obs.matched  = false;
%         end
        % ----- ACTIVE-SEARCH FOR REAL IMAGES IMPLEMENTATION ----- 
        % --------------------------------------------------------
    otherwise
        
        error('??? Unknown Raw data type ''%s''.',Raw.type)
        
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

