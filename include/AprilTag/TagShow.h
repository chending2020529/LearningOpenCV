
#ifndef TAGSHOW_H
#define TAGSHOW_H
#include <iomanip>
#include <cstring>
#include <vector>
#include <list>
#include <chrono>
//#include <pangolin/pangolin.h>
#include <eigen3/Eigen/Core>
#include <cmath>
#include "opencv2/opencv.hpp"
// April tags detector and various families that can be selected by command line option
#include "TagDetector.h"
#include "Tag16h5.h"
#include "Tag25h7.h"
#include "Tag25h9.h"
#include "Tag36h9.h"
#include "Tag36h11.h"
#ifndef PI
const double PI = 3.14159265358979323846;
#endif
const double TWOPI = 2.0*PI;
typedef Eigen::Matrix<double, 6,1> Vector6d;

class TagShow {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        TagShow(cv::Mat & insMatrix,cv::Mat &distCoeff, double tagSize,std::string tagType="36h11");
	~TagShow();

	std::vector<double> getCameraPose(){return mPose;};
	std::vector<Vector6d> getCamMatrix(){return mCameraPoseMatrix;};
	Vector6d localization(AprilTags::TagDetection &detections);
	bool localImg(cv::Mat &img,double tagDis,int tagNum);
	bool processImage(cv::Mat& image);
  bool getRTPos(cv::Mat &image,double xDis,double yDis,double zDis, Vector6d &);
  std::pair<float,float> getTagCenter(){return mCenter;};
  /**
   * Convert rotation matrix to Euler angles
   */
      void wRo_to_euler(const Eigen::Matrix3d& wRo, double& yaw, double& pitch, double& roll);
private:
AprilTags::TagDetector* m_tagDetector;
AprilTags::TagCodes m_tagCodes;
std::vector<Vector6d> mCameraPoseMatrix;
cv::Mat mInsMatrix;
cv::Mat mDistCoeff;
double m_tagSize;
std::vector<double> mPose;
std::pair<float,float> mCenter{0,0};

/**
 * Normalize angle to be within the interval [-pi,pi].
 */
	inline double standardRad(double t) {
  if (t >= 0.) {
    t = fmod(t+PI, TWOPI) - PI;
  } else {
    t = fmod(t-PI, -TWOPI) + PI;
  }
  return t;
}



}; // Demo

#endif
