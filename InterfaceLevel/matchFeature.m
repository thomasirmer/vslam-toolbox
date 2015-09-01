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
            
            hold on
            
            % --------------------------------------------------------
            % normalize observation
            obs_norm = Obs.exp.e / max(Obs.exp.e(1:2));
            
            % --------------------------------------------------------
            % homogenous coordinates of observation
            obs_l = obs_norm(1);
            obs_m = obs_norm(2);
            obs_n = obs_norm(3);
            
            % --------------------------------------------------------
            % calculate neighborhood range (n-sigma);
            n = 20;
            l_sigma = n * sqrt(Obs.exp.E(1,1));
            m_sigma = n * sqrt(Obs.exp.E(2,2));
            n_sigma = n * sqrt(Obs.exp.E(3,3));
            
            % --------------------------------------------------------
            % homogenous coordinates of raw lines
            nLines = size(rawDataLmks.lines, 1);
            for i = 1:nLines
                [rawDataLmks.hmgCoords(:,i), ~] = seg2hmgLin(rawDataLmks.coord(:,i));
            end
            
            % --------------------------------------------------------
            % check all possible lines
            j = 1;
            possibleMatches = zeros(nLines, size(rawDataLmks.lines, 2) + 1);
            for i = 1:nLines
                
                % --------------------------------------------------------
                % normalize raw line
                raw_norm = rawDataLmks.hmgCoords(:,i) / max(rawDataLmks.hmgCoords(1:2,i));
                
                % --------------------------------------------------------
                % homogenous coordinates of raw line
                raw_l = raw_norm(1);
                raw_m = raw_norm(2);
                raw_n = raw_norm(3);
                
                % --------------------------------------------------------
                % normalization factors
                %norm_l = obs_l / raw_l;
                %norm_m = obs_m / raw_m;
                %norm_n = obs_n / raw_n;
                
                % --------------------------------------------------------
                % check if (m, n) in n-sigma range
                %raw_norm_l = raw_l * norm_l;
                %raw_norm_m = raw_m * norm_l;
                %raw_norm_n = raw_n * norm_l;
                
                if (raw_l > obs_l - l_sigma && raw_l < obs_l + l_sigma && ...
                    raw_m > obs_m - m_sigma && raw_m < obs_m + m_sigma && ...
                    raw_n > obs_n - n_sigma && raw_n < obs_n + n_sigma)
                    
                    possibleMatches(j,1:80) = rawDataLmks.lines(i,:);
                    possibleMatches(j,81) = i;
                    plotSegLine(rawDataLmks.coord(:,i), 'blue--*');
                    j = j + 1;
                end
            end
            
            plotHmgLine(Obs.exp.e, 'green');
            plotSegLine(Obs.meas.y, 'magenta--*');
                        
            % --------------------------------------------------------
            % select only possible matches
            possibleMatches = possibleMatches(1:j-1,:);
            
            %matchedLine = MatchLines(Obs.meas.line, possibleMatches);
            MatchLinesByDesc(Obs.meas.line, possibleMatches);
            
            Obs.meas.y   = zeros(size(Obs.meas.y)); % TODO: Get coords of matched line
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

