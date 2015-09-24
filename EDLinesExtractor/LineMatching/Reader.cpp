#include "Reader.h"


Reader::Reader(void)
{
}
void Reader::readCamMat(std::string fileName, cv::Mat &camMat, cv::Mat &distCoeff){
	cv::FileStorage fs(fileName,cv::FileStorage::READ);
	cv::FileNode d = fs["cameraMatrix"];
	cv::read(d,camMat);

	d = fs["distCoeffs"];
	cv::read(d,distCoeff);
	fs.release();
}

void Reader::readPoseFile(const std::string &inputFile, std::vector<cv::Mat> &Rots, std::vector<cv::Mat> &Ts, std::vector<cv::Mat> &Ks, std::vector<std::string> &imageFiles){
	cv::FileStorage fs;
	fs.open(inputFile, cv::FileStorage::READ);
	if (fs.isOpened()){
		int numberPoses;
		fs["numberPoses"] >> numberPoses;
		Rots.resize(numberPoses);
		Ts.resize(numberPoses);
		Ks.resize(numberPoses);
		imageFiles.resize(numberPoses);
		for (int i = 0; i < numberPoses; i++){
			std::stringstream ss1;
			ss1 << "cam_" << i;
			std::stringstream ss2;
			ss2 << "rotation_" << i;
			std::stringstream ss3;
			ss3 << "translation_" << i;
			std::stringstream ss4;
			ss4 << "imagePath_" << i;

			fs[ss1.str()] >> Ks[i];
			fs[ss2.str()] >> Rots[i];
			fs[ss3.str()] >> Ts[i];
			fs[ss4.str()] >> imageFiles[i];
		}
	}
	else{
		imageFiles.clear();
		Ks.clear();
		Ts.clear();
		Rots.clear();
	}
}

void Reader::readTrainConfig(const std::string &inputFile,int &winSizeW, int &winSizeH, int &blockSize,int &blockSizeH, int &blockStride,int &blockStrideH, int &cellSize,int &cellSizeH, int &nBins, 
							 std::string &posSampleFiles,std::string &negSampleFiles,
							 std::string &featuresFile,std::string &svmModelFile,
							 std::string &descriptorVectorFile, 
							 float &svm_t, float &svm_g,float &svm_c,
							 std::string &svmLight_learnExe){
								 cv::FileStorage fs;
								 fs.open(inputFile, cv::FileStorage::READ);


								 fs ["winSizeW"] >>winSizeW;
								 fs ["winSizeH"] >>winSizeH;
								 fs ["blockSize"] >>blockSize;
								 fs ["blockSizeH"] >>blockSizeH;
								 fs ["blockStride"] >>blockStride;
								 fs ["blockStrideH"] >>blockStrideH;
								 fs ["cellSize"] >>cellSize;
								 fs ["cellSizeH"] >>cellSizeH;
								 fs ["nBins"] >>nBins;
								 fs ["posSampleFiles"] >>posSampleFiles;
								 fs ["negSampleFiles"] >>negSampleFiles;
								 fs ["featuresFile"] >>featuresFile;
								 fs ["svmModelFile"] >>svmModelFile;
								 fs ["descriptorVectorFile"] >>descriptorVectorFile;
								 fs ["svm_t"] >>svm_t;
								 fs ["svm_g"] >>svm_g;
								 fs ["svm_c"] >>svm_c;
								 fs ["svmLight_learnExe"] >>svmLight_learnExe;

								 fs.release();
								 if(cellSizeH==0)
									 cellSizeH=cellSize;
								 if(blockSizeH==0)
									 blockSizeH=blockSize;
								 if(blockStrideH==0)
									 blockStrideH=blockStride;
}
void Reader::readTestConfig( std::string &inputFile,int &winSizeW,  int &winSizeH, 
							int &blockSize,int &blockSizeH,  int &blockStride,int &blockStrideH,  int &cellSize,int &cellSizeH,  int &nBins, 
							std::string &testFilesFile, std::string &testResultFilesFile,
							std::string &testResultsFile,  std::string &svmModelFile,
							std::string &descriptorVectorFile){
								cv::FileStorage fs;
								fs.open(inputFile, cv::FileStorage::READ);


								fs ["winSizeW"] >>winSizeW;
								fs ["winSizeH"] >>winSizeH;
								fs ["blockSize"] >>blockSize;
								fs ["blockSizeH"] >>blockSizeH;
								fs ["blockStride"] >>blockStride;
								fs ["blockStrideH"] >>blockStrideH;
								fs ["cellSize"] >>cellSize;
								fs ["cellSizeH"] >>cellSizeH;
								fs ["nBins"] >>nBins;
								fs ["testFilesFile"] >>testFilesFile;
								fs ["testResultFilesFile"]>> testResultFilesFile;
								fs ["testResultsFile"]>> testResultsFile;
								fs ["svmModelFile"] >>svmModelFile;
								fs ["descriptorVectorFile"]>> descriptorVectorFile;

								fs.release();
								if(cellSizeH==0)
									cellSizeH=cellSize;
								if(blockSizeH==0)
									blockSizeH=blockSize;
								if(blockStrideH==0)
									blockStrideH=blockStride;
}
void Reader::readExtractConfig( std::string &inputFile,int &winSizeW,  int &winSizeH, 
							   std::string &posContoursFile, std::string &negContoursFile,
							   std::string &posTargetsFile,  std::string &negTargetsFile){
								   cv::FileStorage fs;
								   fs.open(inputFile, cv::FileStorage::READ);

								   fs ["winSizeW"] >>winSizeW;
								   fs ["winSizeH"] >>winSizeH;
								   fs ["posContoursFile"] >>posContoursFile;
								   fs ["negContoursFile"] >>negContoursFile;
								   fs ["posTargetsFile"] >>posTargetsFile;
								   fs ["negTargetsFile"] >>negTargetsFile;

								   fs.release();

}
void Reader::readFilesNamesFromFile(const std::string &fileListFile, std::vector<std::string> &files){
	std::ifstream myReadFile;
	myReadFile.open(fileListFile);

	std::string line;
	if (myReadFile.is_open()) {
		while (!myReadFile.eof()) {
			std::getline(myReadFile,line);
			if(line.empty()==false)
				files.push_back(line);
		}
		myReadFile.close();
	}
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

void Reader::readIMUFile( std::string &inputFile, int frame_count, std::vector<std::pair<int,cv::Point3d>> &results){
	std::ifstream myReadFile;
	myReadFile.open(inputFile);
	int samplesCount=0;
	int xyzMod=2;
	double x,y,z;
	unsigned int timestamp,start_timestamp;
	std::string line;
	results.resize(0);
	if (myReadFile.is_open()) {
		while (!myReadFile.eof()) {
			std::getline(myReadFile,line);
			if(line.empty()==false&& (line.find("rVec") != std::string::npos)){
				if(xyzMod==2&&line.find("X") != std::string::npos){
					std::vector<std::string> elems;
					split(line,' ',elems);
					x=atof(elems[3].c_str());
					xyzMod=0;
				}
				else if(xyzMod==0&&line.find("Y") != std::string::npos){
					
					std::vector<std::string> elems;
					split(line,' ',elems);
					y=atof(elems[3].c_str());
					xyzMod=(++xyzMod);
				}
				else if(xyzMod==1&&line.find("Z") != std::string::npos){
					
					std::vector<std::string> elems;
					split(line,' ',elems);
					z=atof(elems[3].c_str());
					timestamp=atof(elems[0].c_str());
					if(samplesCount==0)
						start_timestamp=timestamp;
					results.push_back(std::pair<int,cv::Point3d>(timestamp-start_timestamp,cv::Point3d(x,y,z)));
					xyzMod=(++xyzMod);
					//std::cout<<results[samplesCount].first<<std::endl;
					samplesCount++;
		
				}
				
			}
				
		}
		myReadFile.close();
	}
}


void Reader::readRectangleFile(const std::string &fileListFile, std::vector<std::string> &files,std::vector<std::vector<cv::Point2d>> &contours){
	std::vector<std::string> filesTmp;
	readFilesNamesFromFile(fileListFile,filesTmp);

	contours.resize(filesTmp.size());
	files.resize(filesTmp.size());
	std::vector<cv::Point2d> cont(4);

	double scale =1.;

	std::vector<cv::Point2d> extremePoints(4);

	for(int i=0;i<filesTmp.size();i++){

		/*extremePoints[0]=cv::Point2d(std::numeric_limits<float>::max(),std::numeric_limits<float>::max());
		extremePoints[1]=cv::Point2d(std::numeric_limits<float>::min(),std::numeric_limits<float>::max());
		extremePoints[2]=cv::Point2d(std::numeric_limits<float>::min(),std::numeric_limits<float>::min());
		extremePoints[3]=cv::Point2d(std::numeric_limits<float>::max(),std::numeric_limits<float>::min());
		*/
		std::string firstFile = filesTmp[i];
		std::vector<std::string> splits1,splits2;
		split(splits1,firstFile,";",split1::empties_ok);

		for(int j=1;j<5;j++){
			split(splits2,splits1[j],",");
			cont[j-1]=cv::Point2d(atof(splits2[0].c_str())/scale,atof(splits2[1].c_str())/scale);
			/*if(cont[j-1].x<extremePoints[0].x && cont[j-1].y<extremePoints[0].y)
			extremePoints[0]=cont[j-1];
			else if(cont[j-1].x>extremePoints[1].x && cont[j-1].y<extremePoints[1].y)
			extremePoints[1]=cont[j-1];
			else if(cont[j-1].x>extremePoints[2].x && cont[j-1].y>extremePoints[2].y)
			extremePoints[2]=cont[j-1];
			else if(cont[j-1].x<extremePoints[3].x && cont[j-1].y>extremePoints[3].y)
			extremePoints[3]=cont[j-1];*/

		}
		std::sort(cont.begin(), cont.end(), sortContourXAsc);
		std::sort(cont.begin(), cont.end()-2, sortContourYAsc);
		std::sort(cont.begin()+2, cont.end(), sortContourYAsc);

		cv::Point2d tmp = cont[1];
		cont[1] = cont[2];
		cont[2] = tmp;
		tmp = cont[2];
		cont[2] = cont[3];
		cont[3] = tmp;

		contours[i]=(cont);
		files[i]=(splits1[0]);
	}
}


void Reader::readFMatFile(const std::string &fMatFile, 
						  std::vector<std::vector<cv::Point2d>> &pointVector1,
						  std::vector<std::vector<cv::Point2d>> &pointVector2,
						  std::vector<std::pair<std::string,std::string>> &filePairs){
							  std::ifstream rfile;
							  rfile.open(fMatFile);

							  std::string fileString,pred_fileString="";

							  if (rfile.is_open()) {
								  while (!rfile.eof()) {
									  std::getline(rfile,fileString);
									  if(std::string::npos == fileString.find("#")){
										  if(fileString.empty()){
											  std::string firstFile,secondFile;
											  rfile >> firstFile; 
											  rfile >> firstFile; 
											  rfile >> secondFile; 
											  rfile >> secondFile; 
											  bool swap=false;
											  if(firstFile.compare(secondFile)>0){
												  swap=true;
												  std::swap(firstFile,secondFile);
											  }
											  filePairs.push_back(std::pair<std::string,std::string>(firstFile,secondFile));
											  std::getline(rfile,fileString);
											  std::getline(rfile,fileString);
											  int count = atoi(fileString.c_str());
											  std::vector<cv::Point2d> points1(count),points2(count);
											  int x,y;
											  for(int i=0;i<count;i++){
												  rfile >> fileString; 
												  rfile >> fileString;
												  x =atof(fileString.c_str());
												  rfile >> fileString;
												  y =atof(fileString.c_str());
												  if(swap)
													  points2[i]=(cv::Point2d(x,y));
												  else
													  points1[i]=(cv::Point2d(x,y));
												  rfile >> fileString; 
												  rfile >> fileString;
												  x =atof(fileString.c_str());
												  rfile >> fileString;
												  y =atof(fileString.c_str());
												  if(swap)
													  points1[i]=(cv::Point2d(x,y));
												  else
													  points2[i]=(cv::Point2d(x,y));
											  }
											  pointVector1.push_back(points1);
											  pointVector2.push_back(points2);
										  }
									  }
								  }
								  rfile.close();
							  }
}

Reader::~Reader(void)
{
}
