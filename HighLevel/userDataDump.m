% USERDATADUMP  User data for SLAMTBDUMP - Example with data dumped from images.
%   This is a particular case of USERDATA.M. It is intended for
%   demonstration of the SLAM toolbox with points on real images, but
%   skipping the image processing part of the algorithm. 
%   
%   Use this file together with userDataDump.m. In userDataDump, you have
%   at the end the data field ExpOpt.dump, with the following options:
%     > 'write': dumps, for each frame, a file of the robot's control
%       signal and one file for each sensor containing identifiers and
%       extracted pixels from images.
%     > 'read': not dump. Read dumped files instead to simulate the case
%       of working with real images. See below for file naming and formats.
%
%   DUMPED FILE FORMATS: see readProcessedImg.m and readControlSignal.m for
%   details
%
%   See also USERDATA, SLAMTBDUMP, WRITECONTROLSIGNAL, WRITEPROCESSEDIMG.

%   Copyright 2009-2010 Joan Sola @ LAAS-CNRS.
%   Copyright 2013 Joan Sola


% Time variables 
%   - sampling time, first and last frames
Time = struct(...
  'dt',                   .1,...          % sampling time, seconds
  'firstFrame',           35,...           % first frame #
  'lastFrame',            800);           % last frame #

% Simulated world
%   - Simulation landmark sets, playground dimensions
World = struct(...
  'points',           thickCloister(-6,6,-6,6,1,7),... % 3d point landmarks - see THICKCLOISTER. 
  'segments',         []);  % 3D segments - see HOUSE. 
    
% Robot things with their controls
%   - each robot's type and initial configuration, and controls.
%   - motion models (add new model strings if you need more):
%       'constVel'    6D Constant velocity model
%       'odometry'    6D Odometry model
%   - See EULERANGLES for orientations specifications.
% Robot{1} = struct(...                     % ODOMETRY EXAMPLE
%   'id',                 1,...           % robot identifier
%   'name',               'Dala',...      % robot name
%   'type',               'atrv',...      % type of robot
%   'motion',             'odometry',...  % motion model
%   'position',           [0;-5;0],...     % robot position in map
%   'orientationDegrees', [0;0;0],...     % orientation, in degrees, roll pitch yaw.
%   'positionStd',        [0;0;0],...     % position error, std
%   'orientationStd',     [0;0;0],...     % orient. error, std, in degrees
%   'dx',                 [.08;0;0],...     % position increment 8
%   'daDegrees',          [0;0;0.9],...     % angle increment, degrees 9
%   'dxStd',              0.005*[1;1;1],...  % odo linear error std
%   'daStd',              0.05*[1;1;1]);     % odo ang error std, degrees

Robot{1} = struct(...                     % CONSTANT VELOCITY EXAMPLE
  'id',                 3,...           % robot identifier
  'name',               'Dala',...      % robot name
  'type',               'atrv',...      % type of robot
  'motion',             'constVel',...  % motion model
  'position',           [0;-5;0],...    % robot position in map
  'orientationDegrees', [0;0;0],...     % orientation, in degrees, roll pitch yaw.
  'positionStd',        [0;0;0],...     % position error, std
  'orientationStd',     [0;0;0],...     % orient. error, std, degrees
  'velocity',           [.8;0;0],...     % lin. velocity
  'angularVelDegrees',  [0;0;9],...    % ang. velocity, in degrees
  'velStd',             [.1;0;0],...     % lin. vel. error, std
  'angVelStd',          [0;0;1],...     % ang. vel. error, std, degrees
  'dv',                 [0;0;0],...     % veolcity increment
  'dwDegrees',          [0;0;0],...     % ang. vel. increment, degrees
  'dvStd',              [.01;.01;.01],...  % vel perturbation std
  'dwStd',              [.1;.1;.1]);    % ang vel pert. std, degrees



% Sensor things 
%   - each sensor's type and parameters, noise, non-measurable prior.
%   - Sensor types (add new type strings if you need more):
%       'pinHole'   Pin-hole camera
%   - See EULERANGLES for orientations specifications.
Sensor{1} = struct(...
  'id',                 1,...           % sensor identifier
  'name',               'Micropix',...  % sensor name
  'type',               'pinHole',...   % type of sensor
  'robot',              1,...           % robot where it is mounted
  'position',           [0;0;.6],...    % position in robot
  'orientationDegrees', [-90;0;-90],... % orientation in robot, roll pitch yaw
  'positionStd',        [0;0;0],...     % position error std
  'orientationStd',     [0;0;0],...     % orient. error std
  'imageSize',          [640;480],...   % image size
  'pixErrorStd',        1.5,...         % pixel error std
  'intrinsic',          [320;240;320;320],... % intrinsic params [u0 v0 au av]
  'distortion',         [],...          % distortion params
  'frameInMap',         false,...       % add sensor frame in slam map?
  'imGrid',             struct(...      % grid for Active Search
    'numCells',         [8;6],...         % number of H and V grid cells
    'skipOuter',        true));           % skip outer cells for initialization?

% Sensor{2} = struct(...
%   'id',                 2,...           % sensor identifier
%   'name',               'Micropix',...      % sensor name
%   'type',               'pinHole',...   % type of sensor
%   'robot',              1,...           % robot where it is mounted
%   'position',           [0;-0.15;.6],...     % position in robot
%   'orientationDegrees', [-90;0;-90],...      % orientation in robot, roll pitch yaw
%   'positionStd',        [0;0;0],...     % position error std
%   'orientationStd',     [1.5;1.5;1.5],...     % orient. error std
%   'imageSize',          [640;480],...   % image size
%   'pixErrorStd',        1.0,...         % pixel error std
%   'intrinsic',          [320;240;320;320],... % intrinsic params
%   'distortion',         [],...          % distortion params
%   'frameInMap',         false);         % add sensor frame in slam map?


% Estimation options 
Opt = struct(...
  'map',              struct(...    % options for the map
    'numLmks',        100,...           % number of 3d landmarks
    'lmkSize',        6),...            % Size of landmark
  'correct',          struct(...    % options for lmk correction
    'reprojectLmks',  true,...          % reproject lmks after active search?
    'reparametrize',  true,...          % reparametrize lmk?
    'nUpdates',       10,...            % max simultaneus updates
    'MD2th',          9,...             % Threshold on Mahalanobis distance squared
    'linTestIdp',     0.1,...           % threshold on IDP linearity test
    'lines',          struct(...        % options for line corrections
      'innType',      'ortDst',...          % innovation type for lines
      'extPolicy',    false,...             % line extending policy ?
      'extSwitch',    10)),...              % extension policy switch point in pixels
  'init',             struct(...    % Options for initialization
    'nbrInits',       [1 1],...        % number of inits [firstFrame, otherFrames]
    'initType',       'idpPnt',...      % Type of lmk to use for init
    'idpPnt',         struct(...        % inverse-distance prior
      'nonObsMean',   0.01,...              % mean of non obs
      'nonObsStd',    0.5),...              % std of non obs
    'plkLin',         struct(...        % Plucker prior
      'nonObsMean',   [.1;0],...            % mean of non obs
      'nonObsStd',    [.25;1])),...         % std of non obs
  'obs',              struct(...    % Observation options
    'lines',          struct(...        % lines options
      'minLength',    20)));                % minimum segment length
        

% Simulation options
%   - random
SimOpt = struct(...                    
  'random',           struct(...    % random generator options
    'newSeed',        true,...          % select new random seed?
    'fixedSeed',      1,...             % random seed for non-random runs
    'seed',           []),...           % current seed
  'obs',              Opt.obs);     % Observation options



% Figure options  
%   - view, projection, video, ellipses.
%   - figure projections - mapProj:
%       'persp'     Perspective
%       'ortho'     Orthographic
%   - 3D figure views - mapView - see MAPOBSERVER.
%       [a,e,f]     Custom azimuth/elevation/FOV vector. Distance automatic
%       [a,e,f,d]   custom az/el/fov/distance vector.
%   - 3D figure predefined views (edit mapObserver.m to create/edit views):
%       'top'       Top view
%       'side'      Side view
%       'view'      Generic view
%       'normal'    Normal view
%   - objects colors - two options for color specification:
%       'rgbcmykw'  1-char predifined Matlab colors
%       [r g b]     RGB color vector. [0 0 0] is black, [1 1 1] is white.
FigOpt = struct(...
  'renderer',       'opengl',...    % renderer
  'rendPeriod',     1,...          % frames to skip for faster processing
  'createVideo',    false,...       % create video sequences?
  'map',            struct(...      % map figure options
    'size',         [320 240],...       % map figure size
    'lims',        struct(...           % 3D playground limits
      'xMin',            -10,...             
      'xMax',             10,...
      'yMin',            -10,...
      'yMax',             10,...
      'zMin',            -10,...
      'zMax',             10),...
    'proj',         'persp',...         % projection of the 3d figure
    'view',         'view',...          % viewpoint of the 3d figure [30 45 40 20]
    'orbit',        [0 0],...           % AZ and EL orbit angle increments
    'showSimLmk',   true,...            % show simulated landmarks?
    'showEllip',    true,...            % show ellipsoids?
    'colors',       struct(...          % map figure colors
      'border',     [1 1 1],...             %   [r g b]      
      'axes',       [0 0 0],...             % with:
      'bckgnd',     [1 1 1],...             %   [0 0 0] black
      'simLmk',     .3*[1 1 1],...          %   [1 1 1] white
      'defPnt',     struct(...              % euclidean point colors
        'mean',     'b',...                     % mean dot
        'ellip',    [.7 .7 1]),...              % ellipsoid
      'othPnt',     struct(...              % other point colors
        'mean',     'r',...                     % mean dot
        'ellip',    [1 .7 .7]),...              % ellipsoid
      'defLin',     struct(...              % Plucker line colors
        'mean',     [0 .8 0],...                % mean line
        'ellip',    [.6 1 .6]),...              % ellipsoid
      'othLin',     struct(...              % Plucker line colors
        'mean',     [.8 0 0],...                % mean line
        'ellip',    [1 .6 .6]),...              % ellipsoid
      'simu',       'b',...                 %   or 'r', 'b', etc.   
      'est',        'g',...                 % estimated robots and sensors
      'ground',     [.8 .8 .8],...          % simulated robots and sensors
      'label',      [.0 .5 0])),...         % landmark ID labels
  'sensor',         struct(...      % sensor figures options
    'size',         [320 240],...       % sensor figure size
    'showEllip',    false,...           % show ellipses?
    'colors',       struct(...          % Sensor figure colors:
      'border',     .8*[1 1 1],...          %    
      'axes',       [0 0 0],...             % 
      'bckgnd',     [1 1 1],...             %
      'raw',        .3*[1 1 1],...          % 
      'defPnt',     struct(...              % Default point colors
        'updated',  'c',...                     % updated
        'predicted','b'),...                    % predicted
      'othPnt',     struct(...              % other point colors
        'updated',  'r',...                     % updated
        'predicted','m'),...                    % predicted
      'defLin',     struct(...              % Default line colors
        'meas',     'b',...                     % measurement
        'mean',     'g',...                     % mean line
        'ellip',    'y'),...                    % ellipsoid
      'othLin',     struct(...              % other line colors
        'meas',     'b',...                     % measurement
        'mean',     'm',...                     % mean line
        'ellip',    'r'),...                    % ellipsoid
      'label',      [.5 .5 .5])));          %


% Experiment options 
%   - site name, series gathered, estimation run number 
ExpOpt = struct(...
  'root',         '~/SLAM/',...         % root directory
  'site',         'img',...             % Name of the site
  'dataRun',      1,...                 % Run # on this site
  'estimateRun',  1,...                 % slam run for data and site
  'lmkTypes',     Opt.init.initType,... % types of landmarks used
  'sensingType',  'mono',...            % sensing mode
  'mappingType',  'single',...          % mapping mode
  'dump',         'read',...            % 'read', 'write'> read/write data files from simulation
  'procImgName',  'procImg-r%02d-s%02d-i%06d.txt',...   % Name format of processed image
  'imgName',      'img-r%02d-s%02d-i%06d.jpg',...       % Name format of image
  'controlName',  'control-r%02d-i%06d.txt');           % Name format of control file



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

