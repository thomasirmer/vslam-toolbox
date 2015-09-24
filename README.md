# Slamtoolbox with real images

## GIT

The whole project is checked in into a git repository. It contains the original slamtoolbox, the EDLinesDetection library and all needed references. If you want to use the code from the EDLinesDetection you have to compile it for your platform. This can simply be done by calling the `make.m` in the root folder of the project. The only thing that you are required to do beforehand is to install the OpenCV interface which can be installed into MATLAB just by calling the `visionSupportPackages` command on the command window. This directly installs the `mexOpenCV`command that is needed to compile the C/C++ files from this project.

## Extensions

### slamtb.m

The SIMULATION part has been replaced with the following steps:

**Before MAIN LOOP:**

* declare dataset you want to use

<!-- code -->

    directory = './pathToMyDataset/';  
    fileExtension = '.png'; % or whatever you use  
    allFiles = dir([directory, '*', fileExtension]);
	
**Within MAIN LOOP:**

* SIMULATION part:
	* Init RAW structure to use images.
	* Read next image from dataset.
	* Detect lines using FindNLineFeatures function (explained later).
* ESTIMATION part:
	* correctKnownLandmarks.m --> matchFeature.m
	* select only those lines for a possible matching whose endpoints are within a n-sigma range of the endpoints of the current observation.
	* find matching line for current observation using MatchLinesByDescriptor function (explained later).
	
#### FindNLineFeatures.cpp

This function detects lines in the given image. It is realized as a C/C++ [MEX-function](http://de.mathworks.com/help/matlab/matlab_external/introducing-mex-files.html?refresh=true) and uses the EDLinesDetector library which is part of this project. It also uses OpenCV data types and functions. In fact all MATLAB arrays that are given to this function gets converted into a OpenCV MAT object at first.

	lines = FindNLineFeatures(image, N); % about 64 is a reasonably large number for N
	
The given `image` must be a MATLAB gray image. The number of lines `N` is supposed to limit the line search for performance reasons. It will give you the `N` longest lines as they are assumed to be the most significant ones.

The return value is a array of `lines` that have been found in the `image`. It is important to keep the whole line array because this is made such the MatchLinesByDescriptor function can reconstruct the lines from this array.
	
### matchFeature.m
	
***What are the possible matches?***

If we would try to match the current observation with all other lines that have been found the algorithm would take ages to complete. So we only want to match with a subset of lines. This subset should consist of those lines that are within a certrain range of the current observation.

***How do we find them?***

Have a look into the slamToolbox documentation which can be found in `slamToolbox.pdf` within the project folder. This describes how to extend the algorithm with real images. Unfortunately it does this only for point landmarks so we have to find another way to determine the possible matches. Nevertheless we follow the same approach.

The idea is to have some kind of indication which we can use to check if other observations are within a certain range and if so we would like to add them to our possible matches. This indication are the endpoints of the lines. That means that we check the following condition for each line:

*Is the first endpoint of the current line near the first endpoint of our current observation **AND** is the second endpoint of the current line near the second endpoint of our current observation?*

***What does "near" mean?***

Luckily the observation structure provides us this information besides the endpoints. The endpoints are stored in `Obs.par.endp.e` and next to it we have `Obs.par.endp.E` which is the covariance matrix of the endpoint coordinates. The covariance decreases over time when the landmarks become more precise.

The diagonal elements of the matrix are the squared variances of the x- and y-coordinates of the endpoints. This means we can use this to span a rectangle around each endpoint and check if the endpoints of the other lines are within this rectangle.

#### MatchLinesByDescriptor.cpp

Like the FindNLineFeatures function this is also realized as a C/C++ [MEX-function](http://de.mathworks.com/help/matlab/matlab_external/introducing-mex-files.html?refresh=true) and uses the EDLinesDetector library. It also uses the OpenCV data types and functions.

	matchedLine = MatchLinesByDesc(currentObservation, allPossibleMatches);
	
This function gets called within `matchFeature.m`. The `matchFeature.m` gets called for every observation. So the MatchLinesByDesc function is constructed such that it compares each observation with all possible matches.

The `currentObservation` is the whole line for this observation as it has been returned by the FindNLineFeatures function. It can be found in `Obs.meas.line`.

The `allPossibleMatches` are the one we have determined beforehand.

The returned value is the 0-based index of the matched line or -1 if no matching could be found.
			
### customUserDataLin.m

This file contains a few adjustments compared to the userData from the simulation

* Enter the correct number of frames of your dataset into the Time struct.
* Use the constVel model for the Robot.
* Set the robots velocity and angularVelocity to 0 and enter a meaningful std. error for the velocity and angularVelocity.
* Enter the correct image size into the Sensor struct.
* Enter the correct intrinsic parameters into the Sensor struct.

## Visual Studio Project

The MEX-functions and all corresponding C/C++ functions are located within a Visual Studio project. Since we like syntax highlighting, auto completion and other smart features every IDE should have, it is way more convenient to code C/C++ in Visual Studio rather than in MATLAB. The project knows all references an header files however you can't compile them in MATLAB. But if you have installed the OpenCV interface MATLAB does the rest. See the following example on how to compile a MEX-function that contains OpenCV code inside MATLAB:

	mexOpenCV -g 'FindNLineFeatures.cpp' 'EDLineDetector.cpp' 'LineDescriptor.cpp'
	
The `-g` option compiles the code in debug mode. So you can debug this function using Visual Studio.

***Debugging***

If you want to debug your C/C++ code that runs in MATLAB you can just use Visual Studio. Go to *Debug --> Attach to process...* and select the running MATLAB process. Set your breakpoint and wait for it to be hit. Don't worry that your breakpoints are marked as invalid within Visual Studio. They work anyway. 