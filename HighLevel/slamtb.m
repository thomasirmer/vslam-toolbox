% SLAMTB  An EKF-SLAM algorithm with simulator and graphics.
%
%   This script performs multi-robot, multi-sensor, multi-landmark 6DOF
%   EKF-SLAM with simulation and graphics capabilities.
%
%   Please read slamToolbox.pdf in the root directory thoroughly before
%   using this toolbox.
%
%   - Beginners should not modify this file, just edit USERDATA.M and enter
%   and/or modify the data you wish to simulate.
%
%   - More advanced users should be able to create new landmark models, new
%   initialization methods, and possibly extensions to multi-map SLAM. Good
%   luck!
%
%   - Expert users may want to add code for real-data experiments. 
%
%   See also USERDATA, USERDATAPNT, USERDATALIN.
%
%   Also consult slamToolbox.pdf in the root directory.

%   Created and maintained by
%   Copyright 2008, 2009, 2010 Joan Sola @ LAAS-CNRS.
%   Copyright 2011, 2012, 2013 Joan Sola.
%   Programmers (for parts of the toolbox):
%   Copyright David Marquez and Jean-Marie Codol @ LAAS-CNRS
%   Copyright Teresa Vidal-Calleja @ ACFR.
%   See COPYING.TXT for full copyright license.

%% OK we start here

tic;

% clear workspace and declare globals
clear
global Map    

%% I. Specify user-defined options - EDIT USER DATA FILE userData.m

% userData;           % user-defined data. SCRIPT.
% userDataPnt;        % user-defined data for points. SCRIPT.
% userDataLin;        % user-defined data for lines. SCRIPT.
customUserDataLin;


%% II. Initialize all data structures from user-defined data in userData.m
% SLAM data
[Rob,Sen,Raw,Lmk,Obs,Tim]     = createSlamStructures(...
    Robot,...
    Sensor,...      % all user data
    Time,...
    Opt);

% Simulation data
[SimRob,SimSen,SimLmk,SimOpt] = createSimStructures(...
    Robot,...
    Sensor,...      % all user data
    World,...
    SimOpt);

% Graphics handles
[MapFig,SenFig]               = createGraphicsStructures(...
    Rob, Sen, Lmk, Obs,...      % SLAM data
    SimRob, SimSen, SimLmk,...  % Simulator data
    FigOpt);                    % User-defined graphic options


%% III. Initialize data logging
% TODO: Create source and/or destination files and paths for data input and
% logs.
% TODO: do something here to collect data for post-processing or
% plotting. Think about collecting data in files using fopen, fwrite,
% etc., instead of creating large Matlab variables for data logging.

% Clear user data - not needed anymore
clear Robot Sensor World Time   % clear all user data

%% IV. Main loop
for currentFrame = Tim.firstFrame : Tim.lastFrame
    
    % 1. SIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % Simulate robots
%     for rob = [SimRob.rob]
% 
%         % Robot motion
%         SimRob(rob) = simMotion(SimRob(rob),Tim);
%         
%         % Simulate sensor observations
%         for sen = SimRob(rob).sensors
% 
%             % Observe simulated landmarks
%             Raw(sen) = simObservation(SimRob(rob), SimSen(sen), SimLmk, SimOpt) ;
% 
%         end % end process sensors
% 
%     end % end process robots

hFigure = figure;

% coordsXY = load('points.dat');
for rob = [Rob.rob]
    for sen = Rob(rob).sensors
        Raw(sen).type = 'image';
        Raw(sen).data = struct('points',[],'segments',[], 'img',[]);
        Raw(sen).data.segments = struct('coord',[],'app',[]);
        Raw(sen).data.points = struct('coord',[],'app',[]);
        
        % -----------------------------------------------------------------
        % REPLACE THIS WITH C++-LINE-DETECTION-CODE
        % -----------------------------------------------------------------
        % selCols = coordsXY(coordsXY(:,1) == currentFrame,:);
        % [sortedUm1,sortedIdx1] = sort(selCols(:,7),1,'descend');
        %          
        % Raw(sen).data.segments.coord = selCols(sortedIdx1,3:6)';
        % Raw(sen).data.segments.app = selCols(sortedIdx1,2)';
        % -----------------------------------------------------------------
        
        if (currentFrame < Tim.lastFrame)       
            % ----- READ IMAGES -----
            imagePath1 = sprintf('./Datasets/aerial/%03d.png', currentFrame);
            imagePath2 = sprintf('./Datasets/aerial/%03d.png', currentFrame + 1);
            image1 = imread(imagePath1);
            image2 = imread(imagePath2);
            
            % get number of channels
            [~, ~, d1] = size(image1);
            [~, ~, d2] = size(image2);
            
            if (d1 > 1) % color image --> convert to gray
                imageLeft = mat2gray(rgb2gray(image1));
            else % grey image
                imageLeft = mat2gray(image1);
            end
            
            if (d2 > 1) % color image --> convert to gray
                imageRight = mat2gray(rgb2gray(image2));
            else % grey image
                imageRight = mat2gray(image2);
            end
            
            % ----- LINE MATCHING -----
            [linesLeft, linesRight, matching] = EDLinesExtractor(imageLeft, imageRight);
            
            matching = sortrows(matching,5);
            
            Raw(sen).data.segments.coord = [linesLeft(matching(:,4)+1,2)' ; linesLeft(matching(:,4)+1,4)' ; linesLeft(matching(:,4)+1,3)' ; linesLeft(matching(:,4)+1,5)'];
            Raw(sen).data.segments.app = matching(:,4)';

            % ---- BEGIN DEBUG ----            
            % ---- plot coordsXY for test purposes
%             hFigure;
%             lines = coordsXY(coordsXY(:,1) == currentFrame, 3:6);
%             x = [lines(:,1)' ; lines(:,3)'];
%             y = [lines(:,2)' ; lines(:,4)'];
%             plot(x,y);
%             axis([min(coordsXY(:,3)) max(coordsXY(:,5)) min(coordsXY(:,4)) max(coordsXY(:,6))]);
%             pause(0.125);

            hFigure;
            imshow(imageLeft);
            hold on;
                        
            % ---- plot all lines
%             x1 = [linesLeft(:,2)' ; linesLeft(:,4)'];
%             y1 = [linesLeft(:,3)' ; linesLeft(:,5)'];
%             
%             x2 = [linesRight(:,2)' ; linesRight(:,4)'];
%             y2 = [linesRight(:,3)' ; linesRight(:,5)'];
            
            % ---- plot matching lines
            x1 = [linesLeft(matching(:,4)+1,2)' ; linesLeft(matching(:,4)+1,4)'];
            y1 = [linesLeft(matching(:,4)+1,3)' ; linesLeft(matching(:,4)+1,5)'];
            
            x2 = [linesRight(matching(:,5)+1,2)' ; linesRight(matching(:,5)+1,4)'];
            y2 = [linesRight(matching(:,5)+1,3)' ; linesRight(matching(:,5)+1,5)'];
            
            plot(x1, y1, 'r');
            plot(x2, y2, 'b');
            
            clear x1 y1 x2 y2;
            hold off;
            % ----- END DEBUG -----
        end
        
        %[coord, app] = EDLinesExtractor(Raw(sen).data.img);
        % -----------------------------------------------------------------
    end
end


% 2. ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process robots
    for rob = [Rob.rob]

        % Robot motion
        % NOTE: in a regular, non-simulated SLAM, this line is not here and
        % noise just comes from the real world. Here, the estimated robot
        % is noised so that the simulated trajectory can be made perfect
        % and act as a clear reference. The noise is additive to the
        % control input 'u'.
        Rob(rob).con.u = SimRob(rob).con.u + Rob(rob).con.uStd.*randn(size(Rob(rob).con.uStd));
        Rob(rob) = motion(Rob(rob),Tim);
        
        Map.t = Map.t + Tim.dt;
               
        % Process sensor observations
        for sen = Rob(rob).sensors

            % Observe knowm landmarks
            [Rob(rob),Sen(sen),Lmk,Obs(sen,:)] = correctKnownLmks( ...
                Rob(rob),   ...
                Sen(sen),   ...
                Raw(sen),   ...
                Lmk,        ...   
                Obs(sen,:), ...
                Opt) ;

            % Initialize new landmarks
            ninits = Opt.init.nbrInits(1 + (currentFrame ~= Tim.firstFrame));
            for i = 1:ninits
                [Lmk,Obs(sen,:)] = initNewLmk(...
                    Rob(rob),   ...
                    Sen(sen),   ...
                    Raw(sen),   ...
                    Lmk,        ...
                    Obs(sen,:), ...
                    Opt) ;
            end

        end % end process sensors

    end % end process robots

    
    % 3. VISUALIZATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if currentFrame == Tim.firstFrame ...
            || currentFrame == Tim.lastFrame ...
            || mod(currentFrame,FigOpt.rendPeriod) == 0
        
        % Figure of the Map:
        MapFig = drawMapFig(MapFig,  ...
            Rob, Sen, Lmk,  ...
            SimRob, SimSen, ...
            FigOpt);
        
        if FigOpt.createVideo
            makeVideoFrame(MapFig, ...
                sprintf('map-%04d.png',currentFrame), ...
                FigOpt, ExpOpt);
        end
        
        % Figures for all sensors
        for sen = [Sen.sen]
            SenFig(sen) = drawSenFig(SenFig(sen), ...
                Sen(sen), Raw(sen), Obs(sen,:), ...
                FigOpt);
            
            if FigOpt.createVideo
                makeVideoFrame(SenFig(sen), ...
                    sprintf('sen%02d-%04d.png', sen, currentFrame),...
                    FigOpt, ExpOpt);
            end
            
        end

        % Do draw all objects
        drawnow;
    end    

    % 4. DATA LOGGING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: do something here to collect data for post-processing or
    % plotting. Think about collecting data in files using fopen, fwrite,
    % etc., instead of creating large Matlab variables for data logging.
   
end

%% V. Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter post-processing code here

toc;

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

