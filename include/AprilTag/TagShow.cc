#include "TagShow.h"
// default constructor
// #include <robokit/utils/error/error.h>
// #include <robokit/utils/logger/logger.h>

TagShow::TagShow(cv::Mat &insMatrix,cv::Mat &distCoeff, double tagSize, std::string tagType) : // default settings, most can be modified through command line options (see below)
																				m_tagDetector(NULL),
																				m_tagCodes(AprilTags::tagCodes36h11),
																				mInsMatrix(insMatrix),
                                                                                mDistCoeff(distCoeff),
																				m_tagSize(tagSize)
{
	if (tagType == "16h5")
	{
		m_tagCodes = AprilTags::tagCodes16h5;
	}
	else if (tagType == "25h7")
	{
		m_tagCodes = AprilTags::tagCodes25h7;
	}
	else if (tagType == "25h9")
	{
		m_tagCodes = AprilTags::tagCodes25h9;
	}
	else if (tagType == "36h9")
	{
		m_tagCodes = AprilTags::tagCodes36h9;
	}
	else if (tagType == "36h11")
	{
		m_tagCodes = AprilTags::tagCodes36h11;
	}

	m_tagDetector = new AprilTags::TagDetector(m_tagCodes);
	// mCameraPoseMatrix << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
}
TagShow::~TagShow()
{
}

Vector6d TagShow::localization(AprilTags::TagDetection &detection)
{
	vector<cv::Point3f> objPts;
	vector<cv::Point2f> imgPts;

	double s = m_tagSize / 2.0;

	cv::Vec4d tmp[5];

    tmp[0] = cv::Vec4d(-s, -s,0, 1);
    tmp[1] = cv::Vec4d( s, -s,0, 1);
    tmp[2] = cv::Vec4d( s, s,0, 1);
    tmp[3] = cv::Vec4d( -s, s,0, 1);

	objPts.push_back(cv::Point3f(tmp[0][0], tmp[0][1], tmp[0][2]));
	objPts.push_back(cv::Point3f(tmp[1][0], tmp[1][1], tmp[1][2]));
	objPts.push_back(cv::Point3f(tmp[2][0], tmp[2][1], tmp[2][2]));
	objPts.push_back(cv::Point3f(tmp[3][0], tmp[3][1], tmp[3][2]));

	imgPts.push_back(cv::Point2f(detection.p[0].first, detection.p[0].second));
	imgPts.push_back(cv::Point2f(detection.p[1].first, detection.p[1].second));
	imgPts.push_back(cv::Point2f(detection.p[2].first, detection.p[2].second));
	imgPts.push_back(cv::Point2f(detection.p[3].first, detection.p[3].second));

	cv::Mat rvec, tvec;
	Eigen::Vector3d translation;
	Eigen::Vector3d m_translation;
	Eigen::Matrix3d rotation;
	// cv::Matx33d cameraMatrix(
	// 		m_fx, 0, m_px,
	// 		0, m_fy, m_py,
	// 		0, 0, 1);
//	cv::Vec4d distParam(0, 0, 0, 0);

    cv::solvePnP(objPts, imgPts, mInsMatrix, mDistCoeff, rvec, tvec, 0, cv::SOLVEPNP_ITERATIVE);
	//solvePnPRansac(objPts, imgPts, cameraMatrix, distParam, rvec, tvec);
	//cv::solvePnP(objPts, imgPts, cameraMatrix, distParam, rvec, tvec);
	// Eigen::Vector3d vec6d;
     cv::Matx33d r;
     cv::Rodrigues(rvec, r);
     Eigen::Matrix3d wRo;
     wRo << r(0, 0), r(0, 1), r(0, 2), r(1, 0), r(1, 1), r(1, 2), r(2, 0), r(2, 1), r(2, 2);
     double yaw,pitch,roll;
     wRo_to_euler(wRo,yaw,pitch,roll);
	// Eigen::Matrix4d T;
	// T.topLeftCorner(3, 3) = wRo;
	// T.col(3).head(3) << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
	// T.row(3) << 0, 0, 0, 1;
    //  LogInfo("pos = "<<tvec.at<double>(0)<<" "<<tvec.at<double>(1)<<" "<<tvec.at<double>(2));
	Vector6d T;
    T << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2), yaw, pitch, roll;
	return T; //
}

void TagShow::wRo_to_euler(const Eigen::Matrix3d &wRo, double &yaw, double &pitch, double &roll)
{
	double sy = sqrt(wRo(0,0)*wRo(0,0)+wRo(1,0)*wRo(1,0));
	bool singular = sy <1e-6;
	if (!singular)
	{
		yaw = standardRad(atan2(wRo(1, 0), wRo(0, 0)));
		double c = cos(yaw);
		double s = sin(yaw);
		pitch = standardRad(atan2(-wRo(2, 0), wRo(0, 0) * c + wRo(1, 0) * s));
		roll = standardRad(atan2(wRo(0, 2) * s - wRo(1, 2) * c, -wRo(0, 1) * s + wRo(1, 1) * c));
	}
	else
	{
		// LogInfo("singular");
		std::cout << "singular" << std::endl;
		yaw = 0;
		pitch = standardRad(atan2(-wRo(2,0),sy));
		roll = standardRad(atan2(-wRo(1,2),wRo(1,1)));
	}
}
bool TagShow::processImage(cv::Mat &image)
{
	bool detectTag = false;
	cv::Mat image_gray;
	cv::cvtColor(image, image_gray, cv::COLOR_BGR2GRAY);

	vector<AprilTags::TagDetection> FirstDetections = m_tagDetector->extractTags(image_gray);

	mPose.clear();
	int count = 0;
	// double mindis = 99999.0;
	mCameraPoseMatrix.clear();
	// LogInfo("tag num = "<<FirstDetections.size());
	for (int i = 0; i < FirstDetections.size(); ++i)
	{
		auto posVec = localization(FirstDetections[i]);
		FirstDetections[i].draw(image);
        mCenter = FirstDetections[i].cxy;
		mCameraPoseMatrix.push_back(posVec);
        //  double dis = sqrt(posVec(0)*posVec(0)+posVec(1)*posVec(1)+posVec(2)*posVec(2));
		//  if(dis <mindis)
		//  {
		// 	 mCameraPoseMatrix=posVec;
		// 	 mindis=dis;
		//  }
		count++;
		detectTag = true;
	    // LogInfo("Detect tag position "<<posVec(0)<<" "<<posVec(1)<<" "<<posVec(2)<<" "<<posVec(3)<<" "<<posVec(4)<<" "<<posVec(5));
	}
	// if (count > 1)
	// {
	// 	LogInfo("Detect more than one tag");
	// }
	
	return detectTag;
	// optionally send tag information to serial port (e.g. to Arduino)
}

bool TagShow::getRTPos(cv::Mat &image,double xDis,double yDis,double zDis,Vector6d & posv6d)
{
    bool detectTag = false;
	cv::Mat image_gray;
	cv::cvtColor(image, image_gray, cv::COLOR_BGR2GRAY);

	vector<AprilTags::TagDetection> FirstDetections = m_tagDetector->extractTags(image_gray);


	// LogInfo("tag num = "<<FirstDetections.size());
	if(FirstDetections.size()<4)
	  return false;

	vector<cv::Point3f> objPts;
	vector<cv::Point2f> imgPts;
    objPts.push_back(cv::Point3f(0.0,-xDis/2.0,yDis/2.0+zDis));
	objPts.push_back(cv::Point3f(0.0,xDis/2.0,yDis/2.0+zDis));
    objPts.push_back(cv::Point3f(0.0,-xDis/2.0,-yDis/2.0+zDis));
    objPts.push_back(cv::Point3f(0.0,xDis/2.0,-yDis/2.0+zDis));
    
	Eigen::Vector2f pos[4];
	for(int i=0;i<FirstDetections.size();++i)
	{
		if(FirstDetections[i].id>3)
		  return false;
		Eigen::Vector2f point((FirstDetections[i].p[0].first+FirstDetections[i].p[1].first+FirstDetections[i].p[2].first+FirstDetections[i].p[3].first)/4.0,
		(FirstDetections[i].p[0].second+FirstDetections[i].p[1].second+FirstDetections[i].p[2].second+FirstDetections[i].p[3].second)/4.0);
		pos[FirstDetections[i].id] = point;
		FirstDetections[i].draw(image);
        detectTag = true;

	}
	imgPts.push_back(cv::Point2f(pos[0](0),pos[0](1)));
	imgPts.push_back(cv::Point2f(pos[1](0),pos[1](1)));
	imgPts.push_back(cv::Point2f(pos[2](0),pos[2](1)));
	imgPts.push_back(cv::Point2f(pos[3](0),pos[3](1)));
    for(int i=0;i<imgPts.size();++i)
	{
		// LogInfo("imgpoint = "<<pos[i](0)<<" "<<pos[i](1));
	}

	cv::Mat rvec, tvec;
	Eigen::Vector3d translation;
	Eigen::Matrix3d rotation;
	cv::Matx33d r;

//	cv::Vec4d distParam(0, 0, 0, 0);

    cv::solvePnP(objPts, imgPts, mInsMatrix, mDistCoeff, rvec, tvec, 0, cv::SOLVEPNP_ITERATIVE);
	cv::Rodrigues(rvec, r);
	Eigen::Matrix3d wRo;
	double yaw,pitch,roll;
	wRo << r(0, 0), r(0, 1), r(0, 2), r(1, 0), r(1, 1), r(1, 2), r(2, 0), r(2, 1), r(2, 2);
	auto wr = wRo.inverse();
	wRo_to_euler(wr,yaw,pitch,roll);
	translation<<tvec.at<double>(0),tvec.at<double>(1),tvec.at<double>(2);
	auto t = -wr*translation;
	posv6d<<t(0),t(1),t(2),yaw,pitch,roll;
    // LogInfo("pos = "<<t(0)<<" "<<t(1)<<" "<<t(2)<<" "<<roll*180.0/M_PI<<" "<<pitch*180.0/M_PI<<" "<<yaw*180.0/M_PI);



	// Eigen::Vector3d pos[4],calib[4],p,q;
	// for (int i = 0; i < FirstDetections.size(); ++i)
	// {
	// 	auto posVec = localization(FirstDetections[i]);
	// 	FirstDetections[i].draw(image);
	// 	Eigen::Vector3d vec;
	// 	vec<<posVec(0),posVec(1),posVec(2);
    //     pos[FirstDetections[i].id]=vec;
	// 	LogInfo("tag id = "<<FirstDetections[i].id<<" "<<posVec(0)<<" "<<posVec(1)<<" "<<posVec(2));
	// 	detectTag = true;
	// }
	// calib[0]<<0,-xDis/2.0,yDis/2.0+zDis;
	// calib[1]<<0,xDis/2.0,yDis/2.0+zDis;
	// calib[2]<<0,-xDis/2.0,-yDis/2.0+zDis;
	// calib[3]<<0,xDis/2.0,-yDis/2.0+zDis;
	// p = (pos[0]+pos[1]+pos[2]+pos[3])/4;
	// q = (calib[0]+calib[1]+calib[2]+calib[3])/4;
	// for(int i=0;i<4;++i)
	// {
	// 	pos[i]= pos[i]-p;
	// 	calib[i] = calib[i]-q;
	// }
	// Eigen::Matrix3d s = Eigen::Matrix3d::Zero();
	// for(int i=0;i<4;++i)
	// {
	// 	s +=pos[i]*calib[i].transpose();
	// }
	// //SVD
    // Eigen::JacobiSVD<Eigen::Matrix3d> svd(s,Eigen::ComputeFullU | Eigen::ComputeFullV);
	// const Eigen::Matrix3d U = svd.matrixU();
	// const Eigen::Matrix3d V = svd.matrixV();
	// //mirror
	// Eigen::Matrix3d remove_mirror{Eigen::Matrix3d::Identity()};
	// remove_mirror(2,2) = (V*U.transpose()).determinant();

	// auto R = V*remove_mirror*U.transpose();
	// auto t = q - R*p;
	// double yaw,pitch,roll;
	// wRo_to_euler(R,yaw,pitch,roll);


	// LogInfo("pos = "<<t(0)<<" "<<t(1)<<" "<<t(2)<<" "<<roll*180.0/M_PI<<" "<<pitch*180.0/M_PI<<" "<<yaw*180.0/M_PI);
    // posv6d<<t(0),t(1),t(2),yaw,pitch,roll;
    
	return detectTag;
}

bool TagShow::localImg(cv::Mat &Img, double tagDis, int tagNum)
{
	bool detectTag = false;
	cv::Mat image_gray;
	cv::cvtColor(Img, image_gray, cv::COLOR_BGR2GRAY);

	vector<AprilTags::TagDetection> Detections = m_tagDetector->extractTags(image_gray);

	if (Detections.size() == 0)
		return detectTag;

	vector<cv::Point3f> objPts;
	vector<cv::Point2f> imgPts;

	double s = m_tagSize / 2.0;

	cv::Vec4d tmp[5];

	for (int i = 0; i < Detections.size(); ++i)
	{
		if (Detections[i].id >= tagNum)
			break;
		tmp[0] = cv::Vec4d(0, -s + Detections[i].id * tagDis, -s, 1);
		tmp[1] = cv::Vec4d(0, s + Detections[i].id * tagDis, -s, 1);
		tmp[2] = cv::Vec4d(0, s + Detections[i].id * tagDis, s, 1);
		tmp[3] = cv::Vec4d(0, -s + Detections[i].id * tagDis, s, 1);
		objPts.push_back(cv::Point3f(tmp[0][0], tmp[0][1], tmp[0][2]));
		objPts.push_back(cv::Point3f(tmp[1][0], tmp[1][1], tmp[1][2]));
		objPts.push_back(cv::Point3f(tmp[2][0], tmp[2][1], tmp[2][2]));
		objPts.push_back(cv::Point3f(tmp[3][0], tmp[3][1], tmp[3][2]));

		imgPts.push_back(cv::Point2f(Detections[i].p[0].first,
									 Detections[i].p[0].second));
		imgPts.push_back(cv::Point2f(Detections[i].p[1].first,
									 Detections[i].p[1].second));
		imgPts.push_back(cv::Point2f(Detections[i].p[2].first,
									 Detections[i].p[2].second));
		imgPts.push_back(cv::Point2f(Detections[i].p[3].first,
									 Detections[i].p[3].second));
	}
	if (objPts.size() == 0 || imgPts.size() == 0 || objPts.size() != imgPts.size())
		return false;
	// std::cout<<"point size = "<<objPts.size()<<" img size = "<<imgPts.size()<<std::endl;

	cv::Mat rvec, tvec;
	Eigen::Vector3d translation;
	Eigen::Vector3d m_translation;
	Eigen::Matrix3d rotation;
	cv::Vec4d distParam(0, 0, 0, 0);

	cv::solvePnP(objPts, imgPts, mInsMatrix, distParam,
				 rvec, tvec, 0, cv::SOLVEPNP_ITERATIVE);
	cv::Matx33d r;
	cv::Rodrigues(rvec, r);
	Eigen::Matrix3d wRo;
	wRo << r(0, 0), r(0, 1), r(0, 2), r(1, 0),
		r(1, 1), r(1, 2), r(2, 0), r(2, 1), r(2, 2);

	Eigen::Matrix4d T;
	T.topLeftCorner(3, 3) = wRo;
	T.col(3).head(3) << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
	T.row(3) << 0, 0, 0, 1;
	auto T_in = T.inverse();
	auto trans = T_in.col(3).head(3);
	auto rota = T_in.block(0, 0, 3, 3);

	mPose.clear();
	double yaw, pitch, roll;
	wRo_to_euler(rota, yaw, pitch, roll);
	mPose.push_back(trans[0]);
	mPose.push_back(trans[1]);
	mPose.push_back(trans[2]);
	mPose.push_back(roll);
	mPose.push_back(pitch);
	mPose.push_back(yaw);
	return true;
}
// The processing loop where images are retrieved, tags detected,
// and information about detections generated
//void TagShow::loop() {

// Eigen::Vector3d translation;
// Eigen::Matrix3d rotation;
// cv::Mat image_gray;

// int frame = 0;
// double last_t = tic();
// while (true) {

// 	rs2::frameset frames;
// 	try {
// 		frames = m_pipe.wait_for_frames(500);
// 	}
// 	catch (const std::exception& e) {
// 		break;
// 	}
// 	// Try to get a frame of a depth image
// 	//rs2::depth_frame depth = frames.get_depth_frame();
// 	auto color = frames.get_color_frame();

// 	// Get the depth frame's dimensions
// 	float width = color.get_width();
// 	float height = color.get_height();

// 	cv::Mat image(cv::Size(width, height), CV_8UC3, (void*)color.get_data(), cv::Mat::AUTO_STEP);

// 	//cv::remap(image, image, g_mapx, g_mapy, cv::INTER_LINEAR);
// 	processImage(image, image_gray);
// 	translation = mCameraPoseMatrix.col(3).head(3);
// 	rotation = mCameraPoseMatrix.block(0, 0, 3, 3);
// 	double alpha, beta, gamma;
// 	wRo_to_euler(rotation, alpha, beta, gamma);//������������ת�Ƕ�

// 	// exit if any key is pressed
// }
// cv::destroyAllWindows();
//}
