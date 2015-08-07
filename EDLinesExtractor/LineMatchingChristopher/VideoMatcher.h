#pragma once
#include "../PairwiseLineMatching.hpp"
#include "LineFeature.h"
#include <set>
#include "Reader.h"

//#define USE_CV_EDLINES
//#define SHOW_DEBUG
//#define USE_CONCURRENCY
#if defined(USE_CONCURRENCY)
#include "../Concurrency/BlockingQueue.h"
template<typename T>class TaskConsumer{
private:
	BlockingQueue<T> *queue;
	int waitTimeMS;
	std::string name;
public:
	TestConsumer(std::string name, BlockingQueue<T> *queue, int msWait) :name(name), queue(queue), waitTimeMS(msWait){ };
	void consumeTask(){
		while (true){
			T data;
			queue->wait_and_pop(data);
			std::cout << name << " popped: " << data << std::endl;
			boost::this_thread::sleep_for(boost::chrono::milliseconds(200));
		}
	}
	~TestConsumer(){};
};
#endif

#if defined(NDEBUG)
	#include "TimeUtilCV.h"
#endif
#if defined(USE_CV_EDLINES)
	#include "opencv2\line_descriptor.hpp"
#endif 



class VideoMatcher
{
private:
	cv::Mat  prv_image;
	std::string prefix, dir, suffix, file_cnt;
	LineDescriptor lineDesc;
	ScaleLines   prv_keypoints;
	PairwiseLineMatching lineMatch;
	std::set<LineFeature> features;
	int startFrame, endFrame;
	bool mustInitUndistort;
	double factor;
	cv::Mat calib, distCoeffs;
	void loadImage(const std::string &fileName, cv::Mat &image);
	void loadNextImage(const int i, cv::Mat &image, const bool undistort=false, const bool useFisheye = false);
	cv::Mat map1, map2;
	cv::Mat remap(const cv::Mat &image,  const bool useFisheye = false);

public:
	VideoMatcher(const std::string &prefix, const std::string &dir, const std::string &suffix, const int startFrame, const int endFrame, const double &factor=1.0, const int octaves=3, const std::string &file_cnt = "%06d", const cv::Mat &calib = cv::Mat(), const cv::Mat &distCoeffs = cv::Mat() );

	~VideoMatcher();

	void perform(const bool undistort = false, const bool useFisheye=false);
	void featuresToFile(const std::string &file);
	/*static cv::Scalar randomColor(cv::RNG& rng); */
	cv::Scalar randomColor(cv::RNG& rng)
	{
		int icolor = (unsigned)rng;
		return cv::Scalar(icolor & 255, (icolor >> 8) & 255, (icolor >> 16) & 255);
}
#if defined(SHOW_DEBUG)
	/*static */ cv::Mat  drawLines(cv::Mat &imgLinesL, cv::Mat &imgLinesR, std::vector<cv::DMatch> &matchResult, ScaleLines &linesL, ScaleLines &linesR, int firstN=999999);
#endif
#if  defined(SHOW_DEBUG) && defined(USE_CV_EDLINES)
	/*static */ cv::Mat  drawLines(cv::Mat &imgLinesL, cv::Mat &imgLinesR, std::vector<cv::DMatch> &matchResult, std::vector<cv::line_descriptor::KeyLine> &linesL, std::vector<cv::line_descriptor::KeyLine> &linesR);
#endif
	

#if defined(USE_CV_EDLINES)
	/*static */void keyLines2scaleLines(const std::vector<cv::line_descriptor::KeyLine> &keylines, ScaleLines &scaleLines/*, bool computeDescriptors = true, LineDescriptor &lineDescriptor=LineDescriptor()*/,  bool clearScaleLines = false){
		if (clearScaleLines)
			scaleLines.clear();
		//ScaleLines sl(keylines.size(),LinesVec(1,OctaveSingleLine()));
		for each (cv::line_descriptor::KeyLine line in keylines){
			scaleLines.push_back(LinesVec(1));
			//scaleLines.back().push_back(OctaveSingleLine());

			OctaveSingleLine &singleLine = scaleLines.back().back();

			singleLine.endPointX = line.endPointX;
			singleLine.endPointY = line.endPointY;
			singleLine.startPointX = line.startPointX;
			singleLine.startPointY = line.startPointY;
			singleLine.lineLength = line.lineLength;
			singleLine.numOfPixels = line.numOfPixels;
			singleLine.octaveCount =line.octave;
			singleLine.direction = line.angle ;
			singleLine.ePointInOctaveX = line.ePointInOctaveX;
			singleLine.ePointInOctaveY = line.ePointInOctaveY;
			singleLine.sPointInOctaveX = line.sPointInOctaveX;
			singleLine.sPointInOctaveY = line.sPointInOctaveY;
		}
		/*if (computeDescriptors){
			lineDescriptor.ComputeLBD_(scaleLines);
		}*/
	}
#endif
};


// inSegment(): determine if a point is inside a segment
//    Input:  a point P, and a collinear segment S
//    Return: 1 = P is inside S
//            0 = P is  not inside S
int
inSegment(cv::Point3f & P, cv::Point3f &SP0, cv::Point3f &SP1)
{
	if (SP0.x != SP1.x) {    // S is not  vertical
		if (SP0.x <= P.x && P.x <= SP1.x)
			return 1;
		if (SP0.x >= P.x && P.x >= SP1.x)
			return 1;
	}
	else {    // S is vertical, so test y  coordinate
		if (SP0.y <= P.y && P.y <= SP1.y)
			return 1;
		if (SP0.y >= P.y && P.y >= SP1.y)
			return 1;
	}
	return 0;
}
// intersect2D_2Segments(): find the 2D intersection of 2 finite segments
//    Input:  two finite segments S1 and S2
//    Output: *I0 = intersect point (when it exists)
//            *I1 =  endpoint of intersect segment [I0,I1] (when it exists)
//    Return: 0=disjoint (no intersect)
//            1=intersect  in unique point I0
//            2=overlap  in segment from I0 to I1
//#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product  (2D)
float perp(cv::Point3f &u, cv::Point3f & v){
	return u.x*v.y - u.y*v.x;
}


bool intersection(const cv::Point3f&  afirst, const cv::Point3f&  asecond, const cv::Point3f&  bfirst, const cv::Point3f&  bsecond, cv::Point3f& ip)
// http://mathworld.wolfram.com/Line-LineIntersection.html
// in 3d; will also work in 2d if z components are 0
{
	cv::Point3f da = asecond - afirst;
	cv::Point3f db = bsecond - bfirst;
	cv::Point3f dc = bfirst - afirst;

	if (dc.dot(da.cross(db)) != 0.0) // lines are not coplanar
		return false;
	cv::Point3f daXdb = da.cross(db);
	float s = (dc.cross(db)).dot(daXdb) / sqrt(daXdb.x*daXdb.x + daXdb.y*daXdb.y + daXdb.z*daXdb.z);
	if (s >= 0.0 && s <= 1.0)
	{
		ip = afirst + da * s;// cv::Point3f(s, s, s);
		return true;
	}

	return false;
}

// dist3D_Line_to_Line(): get the 3D minimum distance between 2 lines
//    Input:  two 3D lines L1 and L2
//    Return: the shortest distance between L1 and L2
float
dist3D_Line_to_Line(const cv::Point3f&  L1P0, const cv::Point3f&  L1P1, const cv::Point3f&  L2P0, const cv::Point3f&  L2P1, cv::Point3f& ip1, cv::Point3f& ip2)
{
	cv::Point3f   u = L1P1 - L1P0;
	cv::Point3f   v = L2P1 - L2P0;
	cv::Point3f   w = L1P0 - L2P0;
	double    a = u.dot(u);         // always >= 0
	double    b = u.dot(v);
	double    c = v.dot(v);         // always >= 0
	double    d = u.dot(w);
	double    e = v.dot(w);
	double    D = a*c - b*b;        // always >= 0
	double    sc, tc;
	const double SMALL_NUM = 0.00000001;
	// compute the line parameters of the two closest points
	if (D < SMALL_NUM) {          // the lines are almost parallel
		sc = 0.0;
		tc = (b>c ? d / b : e / c);    // use the largest denominator
	}
	else {
		sc = (b*e - c*d) / D;
		tc = (a*e - b*d) / D;
	}
	ip1 = L1P0+sc*u;
	ip2 = L2P0+tc*v;
	// get the difference of the two closest points
	cv::Point3f   dP = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)

	return sqrt(dP.x*dP.x + dP.y*dP.y + dP.z*dP.z);   // return the closest distance
}
void checkEnv(cv::Mat &img, cv::Point3f &p){
	int x_r = roundf(p.x), y_r = roundf(p.y);
	ushort min = std::numeric_limits<ushort>::max();
	int min_x=-1, min_y=-1;
	for (int x = std::max(0, x_r - 2); x < x_r + 2 && x < img.cols; x++)
		for (int y = std::max(0, y_r - 2); y < y_r + 2 && y < img.rows; y++)
			if (img.at<ushort>(y, x) > 0 && img.at<ushort>(y, x) < min){
				min = img.at<ushort>(y, x);
				min_x = x;
				min_y = y;
			}
	if (min_x == -1){
		p.z = 0;
		return;
	}
	p.x = min_x;
	p.y = min_y;
	p.z = min;
}
void rotatePoint(const cv::Mat m_quat, const float lx, const float ly, const float lz, float &gx, float &gy, float &gz)
{
	const float &r = m_quat.at<float>(0);
	const float &x = m_quat.at<float>(1);
	const float &y = m_quat.at<float>(2);
	const float &z = m_quat.at<float>(3);
	const float t2 = r*x;
	const float t3 = r*y;
	const float t4 = r*z;
	const float t5 = -x*x;
	const float t6 = x*y;
	const float t7 = x*z;
	const float t8 = -y*y;
	const float t9 = y*z;
	const float t10 = -z*z;
	gx = 2 * ((t8 + t10)*lx + (t6 - t4)*ly + (t3 + t7)*lz) + lx;
	gy = 2 * ((t4 + t6)*lx + (t5 + t10)*ly + (t9 - t2)*lz) + ly;
	gz = 2 * ((t7 - t3)*lx + (t2 + t9)*ly + (t5 + t8)*lz) + lz;
}
void composePoint(const cv::Mat m_quat, const cv::Mat m_coords, const float lx, const float ly, const float lz, float &gx, float &gy, float &gz)
{
	rotatePoint(m_quat, lz, ly, lz, gx, gy, gz);
	gx += m_coords.at<float>(0);
	gy += m_coords.at<float>(1);
	gz += m_coords.at<float>(2);
}

void Quaternion_To_Rotation(const cv::Mat &quat, cv::Mat &rot)
{
	const float &quat0 = quat.at<float>(0, 0);
	const float &quat1 = quat.at<float>(1, 0);
	const float &quat2 = quat.at<float>(2, 0);
	const float &quat3 = quat.at<float>(3, 0);

	rot.create(3, 3, CV_32FC1);
	float denom = (quat0 * quat0) + (quat1 * quat1) + (quat2 * quat2) + (quat3 * quat3);
	rot.at<float>(0, 0) = ((quat0 * quat0) + (quat1 * quat1) - (quat2 * quat2) - (quat3 * quat3)) / denom;
	rot.at<float>(0, 1) = ((2 * quat1 * quat2) - (2 * quat0 * quat3)) / denom;
	rot.at<float>(0, 2) = ((2 * quat1 * quat3) + (2 * quat0 * quat2)) / denom;

	rot.at<float>(1, 0) = ((2 * quat1 * quat2) + (2 * quat0 * quat3)) / denom;
	rot.at<float>(1, 1) = ((quat0 * quat0) - (quat1 * quat1) + (quat2 * quat2) - (quat3 * quat3)) / denom;
	rot.at<float>(1, 2) = ((2 * quat2 * quat3) - (2 * quat0 * quat1)) / denom;

	rot.at<float>(2, 0) = ((2 * quat1 * quat3) - (2 * quat0 * quat2)) / denom;
	rot.at<float>(2, 1) = ((2 * quat2 * quat3) + (2 * quat0 * quat1)) / denom;
	rot.at<float>(2, 2) = ((quat0 * quat0) - (quat1 * quat1) - (quat2 * quat2) + (quat3 * quat3)) / denom;

}
float HornMethod(
	const std::vector<cv::Point3f>             &A, const std::vector<cv::Point3f>             &B,
	cv::Mat                   &R, cv::Mat &t,                             // The output vector
	bool forceScaleToUnity)
{

	std::vector<cv::Point3f> A_ = A;
	std::vector<cv::Point3f> B_ = B;

	cv::Mat outVector(7, 1, CV_32FC1);

	// Compute the centroids
	cv::Point3f        cL(0, 0, 0), cR(0, 0, 0);

	const size_t nMatches = A.size();
	//  ASSERT_EQUAL_(input.size()%6, 0)

	for (unsigned int i = 0; i < nMatches; i++)
	{
		cL.x += B[i].x;
		cL.y += B[i].y;
		cL.z += B[i].z;

		cR.x += A[i].x;
		cR.y += A[i].y;
		cR.z += A[i].z;
	}

	//ASSERT_ABOVE_(nMatches,0)
	const double F = 1.0 / nMatches;
	cL *= F;
	cR *= F;

	cv::Mat S = cv::Mat::zeros(3, 3, CV_32FC1); // S.zeros(); // Zeroed by default

	// Substract the centroid and compute the S matrix of cross products
	for (unsigned int i = 0; i < nMatches; i++)
	{
		B_[i].x -= cL.x;
		B_[i].y -= cL.y;
		B_[i].z -= cL.z;

		A_[i].x -= cR.x;
		A_[i].y -= cR.y;
		A_[i].z -= cR.z;

		S.at<float>(0, 0) += A_[i].x*B_[i].x;
		S.at<float>(0, 1) += A_[i].x*B_[i].y;
		S.at<float>(0, 2) += A_[i].x*B_[i].z;

		S.at<float>(1, 0) += A_[i].y*B_[i].x;
		S.at<float>(1, 1) += A_[i].y*B_[i].y;
		S.at<float>(1, 2) += A_[i].y*B_[i].z;

		S.at<float>(2, 0) += A_[i].z*B_[i].x;
		S.at<float>(2, 1) += A_[i].z*B_[i].y;
		S.at<float>(2, 2) += A_[i].z*B_[i].z;
	}

	// Construct the N matrix
	cv::Mat N = cv::Mat::zeros(4, 4, CV_32FC1);; // N.zeros(); // Zeroed by default

	N.at<float>(0, 2) = S.at<float>(0, 0) + S.at<float>(1, 1) + S.at<float>(2, 2);
	N.at<float>(0, 0) = S.at<float>(0, 0) + S.at<float>(1, 1) + S.at<float>(2, 2);
	N.at<float>(0, 1) = S.at<float>(1, 2) - S.at<float>(2, 1);
	N.at<float>(0, 2) = S.at<float>(2, 0) - S.at<float>(0, 2);
	N.at<float>(0, 3) = S.at<float>(0, 1) - S.at<float>(1, 0);

	N.at<float>(1, 0) = N.at<float>(0, 1);
	N.at<float>(1, 1) = S.at<float>(0, 0) - S.at<float>(1, 1) - S.at<float>(2, 2);
	N.at<float>(1, 2) = S.at<float>(0, 1) + S.at<float>(1, 0);
	N.at<float>(1, 3) = S.at<float>(2, 0) + S.at<float>(0, 2);

	N.at<float>(2, 0) = N.at<float>(0, 2);
	N.at<float>(2, 1) = N.at<float>(1, 2);
	N.at<float>(2, 2) = -S.at<float>(0, 0) + S.at<float>(1, 1) - S.at<float>(2, 2);
	N.at<float>(2, 3) = S.at<float>(1, 2) + S.at<float>(2, 1);

	N.at<float>(3, 0) = N.at<float>(0, 3);
	N.at<float>(3, 1) = N.at<float>(1, 3);
	N.at<float>(3, 2) = N.at<float>(2, 3);
	N.at<float>(3, 3) = -S.at<float>(0, 0) - S.at<float>(1, 1) + S.at<float>(2, 2);

	// q is the quaternion correspondent to the greatest eigenvector of the N matrix (last column in Z)
	cv::Mat Z, D;
	cv::Mat v;

	cv::eigen(N, Z, D);
	v = Z.col(Z.cols - 1);

	// ASSERTDEB_( fabs( sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3] ) - 1.0 ) < 0.1 );

	// Make q_r > 0
	if (v.at<float>(0) < 0){
		v.at<float>(0) *= -1;
		v.at<float>(1) *= -1;
		v.at<float>(2) *= -1;
		v.at<float>(3) *= -1;
	}

	cv::Mat q(4, 1, CV_32FC1);          // Create a pose rotation with the quaternion
	for (unsigned int i = 0; i < 4; i++)                    // insert the quaternion part
		outVector.at<float>(i + 3) = q.at<float>(i) = v.at<float>(i);

	// Compute scale
	double  num = 0.0;
	double  den = 0.0;
	for (unsigned int i = 0; i < nMatches; i++)
	{
		num += pow(A_[i].x, 2) + pow(A_[i].y, 2) + pow(A_[i].z, 2);
		den += pow(B_[i].x, 2) + pow(B_[i].y, 2) + pow(B_[i].z, 2);
	} // end-for

	// The scale:
	double s = std::sqrt(num / den);

	// Enforce scale to be 1
	if (forceScaleToUnity)
		s = 1.0;


	cv::Mat pp = cv::Mat::zeros(3, 1, CV_32FC1);
	cv::Mat m_coord = cv::Mat::zeros(3, 1, CV_32FC1);
	composePoint(q, m_coord, cL.x, cL.y, cL.z, pp.at<float>(0), pp.at<float>(1), pp.at<float>(2));
	pp *= s;

	outVector.at<float>(0) = cR.x - pp.at<float>(0);     // X
	outVector.at<float>(1) = cR.y - pp.at<float>(1);     // Y
	outVector.at<float>(2) = cR.z - pp.at<float>(2);     // Z
	//std::cout<<"outVector"<<outVector<<std::endl;
	cv::Mat rot;
	Quaternion_To_Rotation(outVector(cv::Rect(0, 3, 1, 4)), rot);
	//std::cout<<"Rot:\n"<<rot<<std::endl;
	//std::cout<<outVector(cv::Rect(0,0,1,3))<<std::endl;
	//std::cout<<(rot*s)*cv::Mat(B[0])+outVector(cv::Rect(0,0,1,3))<<std::endl;
	rot.copyTo(R);
	outVector(cv::Rect(0, 0, 1, 3)).copyTo(t);
	return s;

}

void calcamposeLines(cv::Mat &XXc, cv::Mat &XXw, cv::Mat &R2, cv::Mat &t2){
	const int n = XXc.cols;
	const int n_half = n/2;
	cv::Mat A(3,n,  CV_64FC1);
	cv::Mat B(3,n,  CV_64FC1);
	for (int i = 0; i < n; i+=2){
		const cv::Mat &s = XXc.col(i);
		const cv::Mat &e = XXc.col(i+1);
		cv::Mat v1 = s - e;
		cv::Mat v2 = e - s;
		A(cv::Rect(i, 0, 1, 3)) = v1*(1. / cv::norm(v1));
		A(cv::Rect(i+1, 0, 1, 3)) = v2*(1. / cv::norm(v2));

		const cv::Mat &s2 = XXw.col(i);
		const cv::Mat &e2 = XXw.col(i + 1);
		v1 = s2 - e2;
		v2 = e2 - s2;
		B(cv::Rect(i, 0, 1, 3)) = v1*(1. / cv::norm(v1));
		B(cv::Rect(i + 1, 0, 1, 3)) = v2*(1. / cv::norm(v2));
	}
	cv::Mat E = A*B.t();
	cv::SVD svd(E);
	cv::Mat R = svd.vt.t()*svd.u.t();
}
void calcamposeV(const std::vector<cv::Point3f> &XXc, const std::vector<cv::Point3f> &XXw, cv::Mat &R2, cv::Mat &t2){
	const int n = XXc.size();
	cv::Mat X(XXw, false);//%B
	cv::Mat Y(XXc,false);//%A
	cv::Mat X_ = X.reshape(1, n);// (XXw);//%B
	cv::Mat Y_ = Y.reshape(1, n);// (XXw);//%B
	X_.convertTo(X_, CV_64FC1);
	Y_.convertTo(Y_, CV_64FC1);
	cv::Mat K = cv::Mat::eye(n, n, CV_64FC1) - cv::Mat::ones(n, n, CV_64FC1) / n;

	const int m = X.channels();
	cv::Mat ux = cv::Mat::zeros(m, 1, CV_64FC1);
	cv::Mat uy = cv::Mat::zeros(m, 1, CV_64FC1);
	cv::Scalar ux_ =cv::mean(X);
	for (int i = 0; i < m; i++)
		ux.at<double>(i, 0) = ux_[i];
	cv::Scalar uy_ = cv::mean(Y);
	for (int i = 0; i < m; i++)
		uy.at<double>(i, 0) = uy_[i];

	cv::Mat XK = X_.t()*K;
	cv::Mat ones31 = cv::Mat::ones(1, XK.cols, CV_64FC1);
	for (int i = 0; i<XK.rows; i++)
		for (int j = 0; j<XK.cols; j++)
			XK.at<double>(i, j) *= XK.at<double>(i, j);

	cv::Mat suma = cv::Mat::zeros(1, XK.cols, CV_64FC1);
	for (int i = 0; i<XK.cols; i++)
		suma.at<double>(i) = cv::sum(XK.col(i))[0];

	double sigmx2 = mean(suma)[0];
	cv::Mat SXY = Y_.t()*K*X_ / n;
	cv::SVD svd(SXY);
	cv::Mat S = cv::Mat::eye(3, 3, CV_64FC1);
	if (determinant(SXY)<0)
		S.at<double>(2, 2) = -1;
	R2 = svd.u*S*svd.vt;

	cv::Mat W = cv::Mat::zeros(3, 3, CV_64FC1);
	W.at<double>(0, 0) = svd.w.at<double>(0);
	W.at<double>(1, 1) = svd.w.at<double>(1);
	W.at<double>(2, 2) = svd.w.at<double>(2);
	cv::Mat WS = W*S;
	double c2 = (trace(WS) / sigmx2)[0];
	t2 = uy - c2*R2*ux;

	X = R2.col(0);
	Y = R2.col(1);
	cv::Mat minusZ = -R2.col(2);
	if (norm(X.cross(Y) + minusZ)>2e-2)
		minusZ.copyTo(R2.col(2));
}
void calcampose(cv::Mat &XXc, cv::Mat &XXw, cv::Mat &R2, cv::Mat &t2){
	int n = XXc.cols;
	cv::Mat X(XXw);//%B
	cv::Mat Y(XXc);//%A

	cv::Mat K = cv::Mat::eye(n, n, CV_64FC1) - cv::Mat::ones(n, n, CV_64FC1) / n;

	cv::Mat ux = cv::Mat::zeros(X.rows, 1, CV_64FC1);
	for (int i = 0; i<X.rows; i++)
		ux.at<double>(i) = cv::mean(X.row(i))[0];
	cv::Mat uy = cv::Mat::zeros(Y.rows, 1, CV_64FC1);
	for (int i = 0; i<Y.rows; i++)
		uy.at<double>(i) = cv::mean(Y.row(i))[0];


	cv::Mat XK = X*K;
	cv::Mat ones31 = cv::Mat::ones(1, XK.cols, CV_64FC1);
	for (int i = 0; i<XK.rows; i++)
		for (int j = 0; j<XK.cols; j++)
			XK.at<double>(i, j) *= XK.at<double>(i, j);

	cv::Mat suma = cv::Mat::zeros(1, XK.cols, CV_64FC1);
	for (int i = 0; i<XK.cols; i++)
		suma.at<double>(i) = cv::sum(XK.col(i))[0];

	double sigmx2 = mean(suma)[0];
	cv::Mat SXY = Y*K*(X.t()) / n;
	cv::SVD svd(SXY);
	cv::Mat S = cv::Mat::eye(3, 3, CV_64FC1);
	if (determinant(SXY)<0)
		S.at<double>(2, 2) = -1;
	R2 = svd.u*S*svd.vt;

	cv::Mat W = cv::Mat::zeros(3, 3, CV_64FC1);
	W.at<double>(0, 0) = svd.w.at<double>(0);
	W.at<double>(1, 1) = svd.w.at<double>(1);
	W.at<double>(2, 2) = svd.w.at<double>(2);
	cv::Mat WS = W*S;
	double c2 = (trace(WS) / sigmx2)[0];
	t2 = cv::Mat(uy) - c2*R2*cv::Mat(ux);

	X = R2.col(0);
	Y = R2.col(1);
	cv::Mat minusZ = -R2.col(2);
	if (norm(X.cross(Y) + minusZ)>2e-2)
		minusZ.copyTo(R2.col(2));
}

// intersect3D_SegmentPlane(): find the 3D intersection of a segment and a plane
//    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
//    Output: *I0 = the intersect point (when it exists)
//    Return: 0 = disjoint (no intersection)
//            1 =  intersection in the unique point *I0
//            2 = the  segment lies in the plane
int intersect3D_SegmentPlane(cv::Point3f &S1, ushort depth, cv::Point3f &I){
	I = depth*S1;// S0 + sI * u;                  // compute segment intersect point
	return 1;
}

void reduceNlongestLines(ScaleLines &lines, int n){
	
	std::sort(lines.begin(), lines.end(), [](const LinesVec & a, LinesVec & b) -> bool			{
		return a[0].lineLength > b[0].lineLength;
	});
	lines.erase(lines.begin() + n, lines.end());
}
int main(int argc, char** argv)
{
	/*const int w_ = 512, h_ = 424;
	char file_str[50];
	//cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE(8, cv::Size(20, 20) );
	cv::RNG randomGen;
	ScaleLines linesL;
	LineDescriptor lineDesc(1);
	PairwiseLineMatching lineMatch;
	cv::Mat img_prv;

	cv::Matx33f K = cv::Matx33f::eye();
	K(0, 0) = 365;
	K(1, 1) = 365;
	K(0, 2) = w_/2;
	K(1, 2) = h_/2;
	
#ifdef _DEBUG
	cv::Mat K_ = cv::Mat::eye(3,3,CV_64FC1);
	K_.at<double>(0, 0) = 365;
	K_.at<double>(1, 1) = 365;
	K_.at<double>(0, 2) = 256;
	K_.at<double>(1, 2) = 212;
#endif
	int startI = 31;
	int64 et1 = cv::getTickCount();
	for (int i = startI; i < 499; i++){
		int64 e1 = cv::getTickCount();
		
		sprintf_s(file_str, "C:\\vsfm\\testDataKinect\\depth_%06d.dat", i);
		//std::cout << "Reading " << file_str << std::endl;
		FILE *f;
		fopen_s(&f, file_str, "rb");
		unsigned short* buffer;
		size_t result;
		buffer = (unsigned short*)malloc(sizeof(unsigned short) * w_ * h_);
		result = fread(buffer, sizeof(unsigned short), w_ * h_, f);
		fclose(f);

		cv::Mat image(h_,w_,  CV_16UC1, buffer);

		ScaleLines linesR;
		lineDesc.OctaveKeyLines(image.clone(), linesR);

#ifdef _DEBUG
		cv::Mat  image_8u;
		cv::blur(image, image_8u, cv::Size(5, 5));
		cv::normalize(image_8u, image_8u, 0, 255, cv::NORM_MINMAX);
		image_8u.convertTo(image_8u, CV_8UC1);
		cv::Mat image_8u3;
		cv::cvtColor(image_8u, image_8u3, CV_GRAY2BGR);
		cv::addWeighted(image_8u3, 0.4, cv::Mat(image.size(), CV_8UC3, cv::Scalar(255, 255, 255)), 1 - .4, 0, image_8u3);
		
		cv::line_descriptor::BinaryDescriptor::Params params;
		params.numOfOctave_ = 1;
		params.reductionRatio = 2;

		cv::Ptr<cv::line_descriptor::BinaryDescriptor> bd = cv::line_descriptor::BinaryDescriptor::createBinaryDescriptor(params);
		std::vector<cv::line_descriptor::KeyLine> keylines1;
		bd->detect(image_8u, keylines1);
		VideoMatcher::keyLines2scaleLines(keylines1, linesR);
#endif

		reduceNlongestLines(linesR,20);
		lineDesc.ComputeLBD_(linesR);

		std::vector<cv::DMatch> matchResult;
		if (i != startI){
			lineMatch.LineMatching(linesL, linesR, matchResult);

			for (int ii = 0; ii < matchResult.size(); ii++)			{
				cv::DMatch &match = matchResult[ii];
				match.distance = linesL[match.queryIdx][0].lineLength + linesR[match.trainIdx][0].lineLength;
				match.distance = abs(linesL[match.queryIdx][0].lineLength - linesR[match.trainIdx][0].lineLength)/match.distance;
			}
			std::sort(matchResult.begin(), matchResult.end(), [](const cv::DMatch & a, const cv::DMatch & b) -> bool			{
				return a.distance < b.distance;
			});


			
		
			cv::Mat X(3, matchResult.size() * 2, CV_64FC1);
			cv::Mat Y(3, matchResult.size() * 2, CV_64FC1);
			std::vector<cv::Point3f> X_(matchResult.size() * 2);
			std::vector<cv::Point3f> Y_(matchResult.size() * 2);
			#ifdef _DEBUG
			std::vector<LineCorrespondence2D3D> lineCorrs;
			lineCorrs.reserve(matchResult.size());
			#endif
			std::vector<std::pair<int, float>> match_dist(matchResult.size());
			for (int ii = 0; ii < matchResult.size(); ii++)
			{
				cv::DMatch &match = matchResult[ii];
				match.distance = std::numeric_limits<float>::max();
				OctaveSingleLine &lineLeft = linesL[match.queryIdx][0];
				OctaveSingleLine &lineRight = linesR[match.trainIdx][0];

				cv::Point3f startL3D(lineLeft.startPointX, lineLeft.startPointY, img_prv.at<ushort>(roundf(std::max(0.f, lineLeft.startPointY)), roundf(std::max(0.f, lineLeft.startPointX))));
				checkEnv(img_prv, startL3D);
				if (startL3D.z < 1)
					continue;
				cv::Point3f startR3D(lineRight.startPointX, lineRight.startPointY, image.at<ushort>(roundf(std::max(0.f, lineRight.startPointY)), roundf(std::max(0.f, lineRight.startPointX))));
				checkEnv(image, startR3D);
				if (startR3D.z < 1)
					continue;
				cv::Point3f endL3D(lineLeft.endPointX, lineLeft.endPointY, img_prv.at<ushort>(roundf(std::max(0.f, lineLeft.endPointY)), roundf(std::max(0.f, lineLeft.endPointX))));
				checkEnv(img_prv, endL3D);
				if (endL3D.z < 1)
					continue;
				cv::Point3f endR3D(lineRight.endPointX, lineRight.endPointY, image.at<ushort>(roundf(std::max(0.f, lineRight.endPointY)), roundf(std::max(0.f, lineRight.endPointX))));
				checkEnv(image, endR3D);
				if (endR3D.z < 1)
					continue;
				cv::Point3f result;
	
				X.at<double>(0, ii * 2) = X_[ii * 2].x = (lineLeft.startPointX - K(0, 2)) / K(0, 0)*startL3D.z;
				X.at<double>(1, ii * 2) = X_[ii * 2].y = (lineLeft.startPointY - K(1, 2)) / K(1, 1)*startL3D.z;
				X.at<double>(2, ii * 2) = X_[ii * 2].z = startL3D.z;
				X.at<double>(0, ii * 2 + 1) = X_[ii * 2+1].x = (lineLeft.endPointX - K(0, 2)) / K(0, 0)*endL3D.z;
				X.at<double>(1, ii * 2 + 1) = X_[ii * 2+1].y = (lineLeft.endPointY - K(1, 2)) / K(1, 1)*endL3D.z;
				X.at<double>(2, ii * 2 + 1) = X_[ii * 2+1].z = endL3D.z;
				Y.at<double>(0, ii * 2) = Y_[ii * 2].x = (lineRight.startPointX - K(0, 2)) / K(0, 0)*startR3D.z;
				Y.at<double>(1, ii * 2) = Y_[ii * 2].y = (lineRight.startPointY - K(1, 2)) / K(1, 1)*startR3D.z;
				Y.at<double>(2, ii * 2) = Y_[ii * 2].z = startR3D.z;
				Y.at<double>(0, ii * 2 + 1) = Y_[ii * 2+1].x = (lineRight.endPointX - K(0, 2)) / K(0, 0)*endR3D.z;
				Y.at<double>(1, ii * 2 + 1) = Y_[ii * 2+1].y = (lineRight.endPointY - K(1, 2)) / K(1, 1)*endR3D.z;
				Y.at<double>(2, ii * 2 + 1) = Y_[ii * 2+1].z = endR3D.z;


				float d2linelength=linesL[match.queryIdx][0].lineLength + linesR[match.trainIdx][0].lineLength;
				//match.distance =;
				match_dist[ii].first = ii;
				match_dist[ii].second = (sqrt(pow(X.at<double>(0, ii * 2) - Y.at<double>(0, ii * 2), 2) +
					pow(X.at<double>(1, ii * 2) - Y.at<double>(1, ii * 2), 2) +
					pow(X.at<double>(2, ii * 2) - Y.at<double>(2, ii * 2), 2))
					+ sqrt(pow(X.at<double>(0, ii * 2 + 1) - Y.at<double>(0, ii * 2 + 1), 2) +
					pow(X.at<double>(1, ii * 2 + 1) - Y.at<double>(1, ii * 2 + 1), 2) +
					pow(X.at<double>(2, ii * 2 + 1) - Y.at<double>(2, ii * 2 + 1), 2))) / d2linelength;

				#ifdef _DEBUG
				Line2D_PNL line2D(lineRight.startPointX, lineRight.startPointY, lineRight.endPointX, lineRight.endPointY);
				Line3D_PNL line3D(startL3D.x, startL3D.y, startL3D.z, endL3D.x, endL3D.y, endL3D.z);
				lineCorrs.push_back(LineCorrespondence2D3D(line2D, line3D));
				X_.push_back(startL3D);
				X_.push_back(endL3D);
				Y_.push_back(startR3D);
				Y_.push_back(endR3D);

				#endif
			}
			std::sort(match_dist.begin(), match_dist.end(), [](const std::pair<int, int> & a, const std::pair<int, int> & b) -> bool			{
				return a.second < b.second;
			});


			int topNmatches = std::min(3, (int)matchResult.size());
			//cv::Mat drawn = VideoMatcher::drawLines(img_prv, image, matchResult, linesL, linesR, topNmatches);


			cv::Mat X_top(3,topNmatches * 2,CV_64FC1);
			cv::Mat Y_top(3,topNmatches * 2, CV_64FC1);
			for (int ii = 0; ii < topNmatches; ii++){
				const int iiInXY = match_dist[ii].first;
				X_top.at<double>(0, ii * 2) = X_[iiInXY * 2].x;
				X_top.at<double>(1, ii * 2) = X_[iiInXY * 2].y;
				X_top.at<double>(2, ii * 2) = X_[iiInXY * 2].z;
				X_top.at<double>(0, ii * 2 + 1) = X_[iiInXY * 2 + 1].x;
				X_top.at<double>(1, ii * 2 + 1) = X_[iiInXY * 2 + 1].y;
				X_top.at<double>(2, ii * 2 + 1) = X_[iiInXY * 2 + 1].z;
				Y_top.at<double>(0, ii * 2) = Y_[iiInXY * 2].x;
				Y_top.at<double>(1, ii * 2) = Y_[iiInXY * 2].y;
				Y_top.at<double>(2, ii * 2) = Y_[iiInXY * 2].z;
				Y_top.at<double>(0, ii * 2 + 1) = Y_[iiInXY * 2 + 1].x;
				Y_top.at<double>(1, ii * 2 + 1) = Y_[iiInXY * 2 + 1].y;
				Y_top.at<double>(2, ii * 2 + 1) = Y_[iiInXY * 2 + 1].z;

			}
			cv::Mat R, t;
			calcampose(X_top, Y_top, R, t);

			#ifdef _DEBUG
			PnL pnl;
			cv::Mat K, Rmat,tvec;
			pnl.computeRotAndPos(K_, lineCorrs, Rmat, tvec);

			pnl.computeRotAndPosRANSAC(K_, lineCorrs, Rmat, tvec,50);
			cv::Point3f p1, p2;
			float dist = dist3D_Line_to_Line(X_[0], X_[1], X_[4], X_[5], p1, p2);
			float dist2 = dist3D_Line_to_Line(Y_[0], Y_[1], Y_[4], Y_[5], p1, p2);
			#endif
			//calcamposeLines(X, Y, R, t);
			//cv::Mat R2, t2;
			//HornMethod(X_, Y_, R2, t2,false);
#ifdef _DEBUG
			int64 e2 = cv::getTickCount();
			double time_elapsed = (double)(e2 - e1) / cv::getTickFrequency();
			std::cout << "time2: " << time_elapsed << std::endl;
#endif
		}

		
		linesL = linesR;
		img_prv = image;
	}
	int64 et2 = cv::getTickCount();
	double timet_elapsed = (double)(et2 - et1) / cv::getTickFrequency();

	std::cout << "timet2: " << timet_elapsed << std::endl;
	*/
	_putenv_s("OPENCV_OPENCL_DEVICE", ":GPU:0");
	// allocate host buffer
	/*cv::Mat host1 = cv::Mat(400, 400, cv::DataType<float>::type);
	cv::Mat host2 = cv::Mat(400, 400, cv::DataType<float>::type);
	cv::Mat host3,host4;
	cv::UMat device1, device2, device3, device4;


	// populate host buffer with random data
	cv::randu(host1, cv::Scalar::all(-20), cv::Scalar::all(20));
	cv::randu(host2, cv::Scalar::all(-20), cv::Scalar::all(20));


	// copy host buffer to device buffer
	// execute arithmetic operations



	host1.copyTo(device1);
	host2.copyTo(device2);
	device3 = cv::UMat::zeros(device1.rows, device1.cols, device1.type());
	cv::gemm(device1, device2, 1, device3, 0, device4);

	TimeUtilCV tu;
	host1.copyTo(device1);
	host2.copyTo(device2);
	device3 = cv::UMat::zeros(device1.rows, device1.cols, device1.type());
	cv::gemm(device1, device2, 1, device3, 0, device4);
	device4.copyTo(host3);
	tu.printElapsed();


	TimeUtilCV tu2;
	cv::Mat host_tmp = cv::Mat::zeros(device1.rows, device1.cols, device1.type());
	cv::gemm(host1, host2, 1, host_tmp, 0, host4);
	tu2.printElapsed();
	*/

	std::string calibFile = "c:\\vsfm\\opteka_massa\\chessboard640\\640photo\\calib.yml";
	cv::Mat camMat, distCoeffs;
	Reader::readCamMat(calibFile, camMat, distCoeffs);


	VideoMatcher vm3("", "C:/vsfm/opteka_massa/testVideo3", ".png", 1, 239,0.7, 5, "%06d", camMat, distCoeffs);
	vm3.perform(true,true);
	vm3.featuresToFile("C:/Users/kropp/Downloads/slamToolbox2/slamToolbox_14_12_05/points.dat");

		/*VideoMatcher vm3("", "C:/vsfm/opteka_massa/testVideo2", ".png", 1, 150,0.7, 3, "%06d", camMat, distCoeffs);
	vm3.perform(true,true);
	vm3.featuresToFile("C:/Users/kropp/Downloads/slamToolbox2/slamToolbox_14_12_05/points.dat");
	*/
	VideoMatcher vm2("rawoutput", "C:/Users/kropp/Downloads/EKF_monoSLAM_1pRANSAC/sequences/ic", ".png", 90, 410, 1.0,5, "%04d");
	vm2.perform();
	vm2.featuresToFile("C:/Users/kropp/Downloads/slamToolbox2/slamToolbox_14_12_05/points.dat");
	
	VideoMatcher vm("","C:\\vsfm\\opteka_massa\\testVideo", ".jpg",10,599,1.0,5);
	vm.perform();
	vm.featuresToFile("C:/Users/kropp/Downloads/slamToolbox2/slamToolbox_14_12_05/points.dat");
}