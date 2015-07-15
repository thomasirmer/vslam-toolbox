#pragma once
#include <vector>
#include <opencv2\opencv.hpp>
class LineFeature
{

public:
	LineFeature();
	~LineFeature();
	int startFrame, endFrame;

	int lastTrainIndex;
	std::vector<std::pair<cv::Point2f, cv::Point2f>> points;
	bool operator<(const LineFeature& rhs) const	{
		if (endFrame == rhs.endFrame)
			return lastTrainIndex<rhs.lastTrainIndex;
		return endFrame<rhs.endFrame;;
	}

	bool operator==(const LineFeature& rhs) const	{
		return (endFrame == rhs.endFrame &&lastTrainIndex == rhs.lastTrainIndex);
	}
};