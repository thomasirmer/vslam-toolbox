#pragma once


#include <fstream>
#include <stdexcept>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "Split.h"

class Reader
{
public:
	inline static bool sortContourXAsc(const cv::Point2d &d1, const  cv::Point2d &d2){
		try{
			return d1.x<d2.x;
		}
		catch(std::exception &e){
			std::cout<<"sortMatchesAscQuery: "<<e.what()<<std::endl;
		}
	}
	inline static bool sortContourXDsc(const cv::Point2d &d1, const  cv::Point2d &d2){
		try{
			return d1.x>d2.x;
		}
		catch(std::exception &e){
			std::cout<<"sortMatchesAscQuery: "<<e.what()<<std::endl;
		}
	}
		inline static bool sortContourYAsc(const cv::Point2d &d1, const  cv::Point2d &d2){
		try{
			return d1.y<d2.y;
		}
		catch(std::exception &e){
			std::cout<<"sortMatchesAscQuery: "<<e.what()<<std::endl;
		}
	}
	Reader(void);
	static void readCamMat(std::string fileName, cv::Mat &camMat, cv::Mat &distCoeff);

	static void readTrainConfig(const std::string &inputFile,int &winSizeW, int &winSizeH, int &blockSize,int &blockSizeH, int &blockStride,int &blockStrideH, int &cellSize,int &cellSizeH, int &nBins, 
							 std::string &posSampleFiles,std::string &negSampleFiles,
							 std::string &featuresFile,std::string &svmModelFile,
							 std::string &descriptorVectorFile, 
							 float &svm_t, float &svm_g,float &svm_c,
							 std::string &svmLight_learnExe );
	static void readTestConfig(std::string &inputFile,int &winSizeW,  int &winSizeH, 
							int &blockSize,int &blockSizeH,  int &blockStride,int &blockStrideH,  int &cellSize,int &cellSizeH,  int &nBins, 
							std::string &testFilesFile, std::string &testResultFilesFile,
							std::string &testResultsFile,  std::string &svmModelFile,
							std::string &descriptorVectorFile);
	static void readExtractConfig( std::string &inputFile,int &winSizeW,  int &winSizeH, 
							 std::string &posContoursFile, std::string &negContoursFile,
							 std::string &posTargetsFile,  std::string &negTargetsFile);
	static void readIMUFile( std::string &inputFile, int frame_count, std::vector<std::pair<int,cv::Point3d>> &results);
	static void readFMatFile(const std::string &fMatFile, 
		std::vector<std::vector<cv::Point2d>> &pointVector1,
		std::vector<std::vector<cv::Point2d>> &pointVector2,
		std::vector<std::pair<std::string,std::string>> &filePairs);
	static void readFilesNamesFromFile(const std::string &inputFile, std::vector<std::string> &files);
	static void readRectangleFile(const std::string &fileListFile, std::vector<std::string> &files,std::vector<std::vector<cv::Point2d>> &contours);

	static void readPoseFile(const std::string &inputFile, std::vector<cv::Mat> &Rots, std::vector<cv::Mat> &Ts, std::vector<cv::Mat> &Ks, std::vector<std::string> &imageFiles);
	~Reader(void);
};

