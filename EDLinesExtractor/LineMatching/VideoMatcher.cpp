#include "VideoMatcher.h"
#include <fstream>

VideoMatcher::VideoMatcher(const std::string &prefix, const std::string &dir, const std::string &suffix, const int startFrame, const int endFrame, const double &factor, const int octaves, const std::string &file_cnt, const cv::Mat &calib, const cv::Mat &distCoeffs)
	: prefix(prefix), dir(dir), file_cnt(file_cnt), suffix(suffix), startFrame(startFrame), endFrame(endFrame), calib(calib), distCoeffs(distCoeffs), mustInitUndistort(true),
	factor(factor)
{


}

/*VideoMatcher::VideoMatcher(const std::string &prefix, const std::string &dir, const std::string &suffix, const int startFrame, const int endFrame, const double &factor, const int octaves, const std::string &file_cnt, const cv::Mat &calib, const cv::Mat &distCoeffs)
: prefix(prefix), dir(dir), file_cnt(file_cnt), suffix(suffix), startFrame(startFrame), endFrame(endFrame), calib(calib), distCoeffs(distCoeffs), mustInitUndistort(true), 
lineDesc(octaves), factor(factor)
{

}*/


void VideoMatcher::loadNextImage(const int i, cv::Mat &image, const bool undistort, const bool useFisheye){
	char base_name[256]; 
	sprintf_s(base_name, (dir+"/" + prefix + file_cnt+suffix).c_str(), i);
	if (undistort)
		image = remap(cv::imread(base_name, CV_LOAD_IMAGE_GRAYSCALE), useFisheye);
	else
		image = cv::imread(base_name, CV_LOAD_IMAGE_GRAYSCALE);
}

#if defined(SHOW_DEBUG)
cv::Mat VideoMatcher::drawLines(cv::Mat &imgLinesL, cv::Mat &imgLinesR, std::vector<cv::DMatch> &matchResult, ScaleLines &linesL, ScaleLines &linesR, int firstN){

	cv::Mat mm(imgLinesL.size().height, imgLinesL.size().width * 2, CV_8UC1);

	imgLinesL.copyTo(mm(cv::Rect(0, 0, imgLinesL.size().width, imgLinesL.size().height)));
	imgLinesR.copyTo(mm(cv::Rect(imgLinesL.size().width, 0, imgLinesR.size().width, imgLinesR.size().height)));

	cv::Mat mmC;

	cv::RNG randomGen(0xFFFFFFFF);
	cv::cvtColor(mm, mmC, cv::COLOR_GRAY2BGR);

	cv::Mat a(mmC.size(), mmC.type());
	a.setTo(cv::Scalar(255, 255, 255));
	double opacity = 0.4;
	cv::addWeighted(mmC, opacity, a, 1 - opacity, 0, mmC);


	for (int i = 0; i < matchResult.size() && i<firstN; i++)
	{
		cv::DMatch &match = matchResult[i];
		cv::Scalar color = randomColor(randomGen);
		cv::line(mmC, cv::Point2f(linesL[match.queryIdx][0].startPointX, linesL[match.queryIdx][0].startPointY), cv::Point2f(linesL[match.queryIdx][0].endPointX, linesL[match.queryIdx][0].endPointY), color, 2);
		cv::line(mmC, cv::Point2f(linesR[match.trainIdx][0].startPointX + imgLinesL.size().width, linesR[match.trainIdx][0].startPointY), cv::Point2f(linesR[match.trainIdx][0].endPointX + imgLinesL.size().width, linesR[match.trainIdx][0].endPointY), color, 2);
		cv::line(mmC, cv::Point2f(linesL[match.queryIdx][0].startPointX, linesL[match.queryIdx][0].startPointY), cv::Point2f(linesR[match.trainIdx][0].startPointX + imgLinesL.size().width, linesR[match.trainIdx][0].startPointY), color, 1);
	}
	return mmC;
}
#endif

#if defined(SHOW_DEBUG)&&defined(USE_CV_EDLINES)
cv::Mat  VideoMatcher::drawLines(cv::Mat &imgLinesL, cv::Mat &imgLinesR, std::vector<cv::DMatch> &matchResult, std::vector<cv::line_descriptor::KeyLine> &linesL, std::vector<cv::line_descriptor::KeyLine> &linesR){

	cv::Mat mm(imgLinesL.size().height, imgLinesL.size().width * 2, CV_8UC1);

	imgLinesL.copyTo(mm(cv::Rect(0, 0, imgLinesL.size().width, imgLinesL.size().height)));
	imgLinesR.copyTo(mm(cv::Rect(imgLinesL.size().width, 0, imgLinesR.size().width, imgLinesR.size().height)));

	cv::Mat mmC;

	cv::RNG randomGen(0xFFFFFFFF);
	cv::cvtColor(mm, mmC, cv::COLOR_GRAY2BGR);

	cv::Mat a(mmC.size(), mmC.type());
	a.setTo(cv::Scalar(255, 255, 255));
	double opacity = 0.4;
	cv::addWeighted(mmC, opacity, a, 1 - opacity, 0, mmC);

	for each (cv::DMatch match in matchResult)
	{
		cv::Scalar color = randomColor(randomGen);
		cv::line(mmC, cv::Point2f(linesL[match.queryIdx].startPointX, linesL[match.queryIdx].startPointY), cv::Point2f(linesL[match.queryIdx].endPointX, linesL[match.queryIdx].endPointY), color, 2);
		cv::line(mmC, cv::Point2f(linesR[match.trainIdx].startPointX + imgLinesL.size().width, linesR[match.trainIdx].startPointY), cv::Point2f(linesR[match.trainIdx].endPointX + imgLinesL.size().width, linesR[match.trainIdx].endPointY), color, 2);
		cv::line(mmC, cv::Point2f(linesL[match.queryIdx].startPointX, linesL[match.queryIdx].startPointY), cv::Point2f(linesR[match.trainIdx].startPointX + imgLinesL.size().width, linesR[match.trainIdx].startPointY), color, 1);
	}
	return mmC;
}
#endif

cv::Mat VideoMatcher::remap(const cv::Mat &image, const bool useFisheye) {

	cv::Mat undistorted;

	cv::Mat newCamMat = calib.clone();
	newCamMat.at<double>(0, 0) *= factor;
	newCamMat.at<double>(1, 1) *= factor;

	if (useFisheye){
		printf("unimplemented code --> useFisheye == true\n");

		//cv::Mat P = cv::Mat::eye(3, 3, CV_32F);

		//cv::fisheye::initUndistortRectifyMap(
		//	calib,  // computed camera matrix
		//	distCoeffs,    // computed distortion matrix
		//	P,     // optional rectification (none) 
		//	newCamMat,//,     //// camera matrix to generate undistorted
		//	cv::Size(image.cols, image.rows),
		//	// size of undistorted
		//	CV_32FC1,      // type of output map
		//	map1, map2);   // the x and y mapping functions
	} else {
		cv::initUndistortRectifyMap(
			calib,  // computed camera matrix
			distCoeffs,    // computed distortion matrix
			cv::Mat(),     // optional rectification (none) 
			cv::Mat(),//newCamMat,//,     //// camera matrix to generate undistorted
			cv::Size(image.cols, image.rows),
			//            image.size(),  // size of undistorted
			CV_32FC1,      // type of output map
			map1, map2);   // the x and y mapping functions

		mustInitUndistort = false;
		calib = newCamMat;
	}
	// Apply mapping functions
	cv::remap(image, undistorted, map1, map2,
		cv::INTER_LINEAR); // interpolation type

	return undistorted;
}

void VideoMatcher::perform(const bool undistort, const bool useFisheye){

	cv::Mat leftImageOriginal, rightImageOriginal;
	std::vector<cv::DMatch> matchResult;

#ifdef USE_CV_EDLINES
	cv::line_descriptor::BinaryDescriptor::Params params;
	params.numOfOctave_ = 1;
	params.reductionRatio = 2;
	cv::Mat descr1, descr2;
	cv::Ptr<cv::line_descriptor::BinaryDescriptor> bd = cv::line_descriptor::BinaryDescriptor::createBinaryDescriptor(params);
	std::vector<cv::line_descriptor::KeyLine> keylines1, keylines2;
#else
	ScaleLines   linesInLeft, linesInRight;
#endif
	
	for (int i = startFrame; i <= endFrame; i++) {
		std::cout << "Performing " << i << "th image..." << std::endl;

		matchResult.clear();
		if (i == startFrame){
			loadNextImage(i, leftImageOriginal, undistort, useFisheye);

#ifdef USE_CV_EDLINES
			bd->detect(leftImage, keylines1);
			bd->compute(leftImageOriginal, keylines1, descr1);
#else
			lineDesc.GetLineDescriptor(leftImageOriginal.clone(), linesInLeft);
#endif
		}
		else{

#ifdef USE_CV_EDLINES
TODO: set prv_keypoints!
#else
			leftImageOriginal = prv_image;
			linesInLeft = prv_keypoints;
#endif
		}

		loadNextImage(i + 1, rightImageOriginal, undistort, useFisheye);

#ifdef USE_CV_EDLINES
		bd->detect(rightImage, keylines2);
		bd->compute(rightImageOriginal, keylines2, descr2);
		/* create a BinaryDescriptorMatcher object */
		cv::Ptr<cv::line_descriptor::BinaryDescriptorMatcher> bdm = cv::line_descriptor::BinaryDescriptorMatcher::createBinaryDescriptorMatcher();

		/* require match */
		bdm->match(descr1, descr2, matchResult);
		/* plot matches */
		cv::Mat outImg;
		std::vector<char> mask(matchResult.size(), 1);
		cv::line_descriptor::drawLineMatches(leftImageOriginal, keylines1, rightImageOriginal, keylines2, matchResult, outImg, cv::Scalar::all(-1), cv::Scalar::all(-1), mask,
			cv::line_descriptor::DrawLinesMatchesFlags::DEFAULT);

#else
		cv::Mat rightImage = rightImageOriginal.clone();
		lineDesc.GetLineDescriptor(rightImage, linesInRight);
		lineMatch.LineMatching(linesInLeft, linesInRight, matchResult);
#endif

#if defined(SHOW_DEBUG) && !defined(USE_CV_EDLINES)
		cv::Mat dl =drawLines(leftImageOriginal, rightImageOriginal, matchResult, linesInLeft, linesInRight);
#elif defined(SHOW_DEBUG) && defined(USE_CV_EDLINES)
		drawLines(leftImageOriginal, rightImageOriginal, matchResult, keylines1, keylines2);
#endif

		std::sort(matchResult.begin(), matchResult.end(), [](cv::DMatch const&d1, cv::DMatch const&d2){
			if (d1.queryIdx == d2.queryIdx)
				return d1.trainIdx < d2.trainIdx;
			else
				return d1.queryIdx < d2.queryIdx; });


			cv::Point2f sP;
			cv::Point2f eP;
			if (i == startFrame){
				for (int m_i = 0; m_i < matchResult.size(); m_i++){
					const cv::DMatch &match_i = matchResult[m_i];
					LineFeature feature;
					feature.startFrame = i;
					feature.endFrame = i + 1;
					feature.lastTrainIndex = match_i.trainIdx;

#if defined(USE_CV_EDLINES)
					feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(keylines1[match_i.queryIdx].startPointX, keylines1[match_i.queryIdx].startPointY), cv::Point2f(keylines1[match_i.queryIdx].endPointX, keylines1[match_i.queryIdx].endPointY)));
					feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(keylines2[match_i.trainIdx].startPointX, keylines2[match_i.trainIdx].startPointY), cv::Point2f(keylines2[match_i.trainIdx].endPointX, keylines2[match_i.trainIdx].endPointY)));
#else
					feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(linesInLeft[match_i.queryIdx][0].startPointX, linesInLeft[match_i.queryIdx][0].startPointY), cv::Point2f(linesInLeft[match_i.queryIdx][0].endPointX, linesInLeft[match_i.queryIdx][0].endPointY)));
					feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(linesInRight[match_i.trainIdx][0].startPointX, linesInRight[match_i.trainIdx][0].startPointY), cv::Point2f(linesInRight[match_i.trainIdx][0].endPointX, linesInRight[match_i.trainIdx][0].endPointY)));

#endif
					features.insert(feature);
				}
			} else {
				for (int m_i = 0; m_i < matchResult.size(); m_i++){
					const cv::DMatch &match_i = matchResult[m_i];

					LineFeature feature;
					feature.startFrame = -1;
					feature.endFrame = i;
					feature.lastTrainIndex = match_i.queryIdx;

					std::set<LineFeature>::iterator a = features.find(feature);
					feature.endFrame = i + 1;
					if (a == features.end()){

#if defined(USE_CV_EDLINES)
						feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(keylines1[match_i.queryIdx].startPointX, keylines1[match_i.queryIdx].startPointY), cv::Point2f(keylines1[match_i.queryIdx].endPointX, keylines1[match_i.queryIdx].endPointY)));
						feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(keylines2[match_i.trainIdx].startPointX, keylines2[match_i.trainIdx].startPointY), cv::Point2f(keylines2[match_i.trainIdx].endPointX, keylines2[match_i.trainIdx].endPointY)));
#else
						feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(linesInLeft[match_i.queryIdx][0].startPointX, linesInLeft[match_i.queryIdx][0].startPointY), cv::Point2f(linesInLeft[match_i.queryIdx][0].endPointX, linesInLeft[match_i.queryIdx][0].endPointY)));
						feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(linesInRight[match_i.trainIdx][0].startPointX, linesInRight[match_i.trainIdx][0].startPointY), cv::Point2f(linesInRight[match_i.trainIdx][0].endPointX, linesInRight[match_i.trainIdx][0].endPointY)));
#endif
						feature.startFrame = i - 1;
					} else {
						feature.startFrame = a->startFrame;
						feature.lastTrainIndex = a->lastTrainIndex;
						feature.points.clear();
						feature.points = a->points;
#if defined(USE_CV_EDLINES)
						feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(keylines2[match_i.trainIdx].startPointX, keylines2[match_i.trainIdx].startPointY), cv::Point2f(keylines2[match_i.trainIdx].endPointX, keylines2[match_i.trainIdx].endPointY)));
#else
						feature.points.push_back(std::pair<cv::Point2f, cv::Point2f>(cv::Point2f(linesInRight[match_i.trainIdx][0].startPointX, linesInRight[match_i.trainIdx][0].startPointY), cv::Point2f(linesInRight[match_i.trainIdx][0].endPointX, linesInRight[match_i.trainIdx][0].endPointY)));
#endif
						features.erase(a);
					}

					feature.lastTrainIndex = match_i.trainIdx;

					features.insert(feature);
				}
			}
#ifdef USE_CV_EDLINES
TODO : set prv_keypoints!
#else
			prv_keypoints = linesInRight;
			prv_image = rightImageOriginal;
#endif
	}

}

void VideoMatcher::featuresToFile(const std::string &file){
	std::ofstream myfile;
	myfile.open(file);

	std::vector<LineFeature> lf_vector;
	lf_vector.reserve(features.size());
	for (std::set<LineFeature>::iterator feature = features.begin(); feature != features.end(); feature++)
		if (feature->endFrame - feature->startFrame >= 4)
			lf_vector.push_back(*feature);

	std::sort(lf_vector.begin(), lf_vector.end(), [](LineFeature const&d1, LineFeature const&d2){
		if (d1.startFrame == d2.startFrame)
			if (d1.endFrame == d2.endFrame)
				return abs(sqrtf(powf(d1.points[0].first.x - d1.points[0].second.x, 2)) + sqrtf(powf(d1.points[0].first.y - d1.points[0].second.y, 2)))
				> abs(sqrtf(powf(d2.points[0].first.x - d2.points[0].second.x, 2)) + sqrtf(powf(d2.points[0].first.y - d2.points[0].second.y, 2)));
			else
				return d1.endFrame > d2.endFrame;
		else
			return d1.startFrame < d2.startFrame; });


		for (int i = startFrame; i < endFrame; i++) {

			for (int j = 0; j < lf_vector.size(); j++){
				LineFeature feature = lf_vector[j];
				if (i >= feature.startFrame
					&&i < feature.endFrame){
						int pointInFeature = i - feature.startFrame;
						const std::pair<cv::Point2f, cv::Point2f> &lf = feature.points[pointInFeature];
						myfile << i - startFrame + 1 << " " << j + 1 << " " << lf.first.x << " " << lf.first.y << " " << lf.second.x << " " << lf.second.y << " " << feature.endFrame - i << "\n"; std::cout << std::endl;
				}



			}
		}
		/*Üstd::vector<std::vector<std::pair<cv::Point2f, cv::Point2f>>> featuresAll;
		std::vector<int> startFrames;
		std::vector<int> endFrames;
		for (int i = startFrame; i<endFrame; i++) {
		for (std::set<LineFeature>::iterator feature = features.begin(); feature != features.end(); feature++){
		if (i >= feature->startFrame
		&&i<feature->endFrame
		&&feature->endFrame - feature->startFrame >= 2){
		featuresAll.push_back(feature->points);
		startFrames.push_back(feature->startFrame);
		endFrames.push_back(feature->endFrame);
		}
		}
		}
		for (int i = startFrame; i<=endFrame; i++) {

		//myfile << "#^BEGIN OF " << i - startFrame << "\n";

		for (int j = 0; j<featuresAll.size(); j++){
		if (i >= startFrames[j] && i <= endFrames[j]){
		const std::pair<cv::Point2f, cv::Point2f> &lf = featuresAll[j][i - startFrames[j]];
		myfile << i - startFrame + 1 << " " << j + 1 << " " << lf.first.x << " " << lf.first.y << " " << lf.second.x << " " << lf.second.y << "\n";
		}
		}
		//myfile << "#!END OF " << i - startFrame << "\n";
		}*/
		myfile.close();

}

VideoMatcher::~VideoMatcher()
{
}