#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <string>
#include "TagShow.h"
#include <iostream>
// #include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/Core>
// #include <opencv2/core/eigen.hpp>

using namespace std;
using namespace cv;

class AgvCalibration
{
public:
    AgvCalibration(double atsize, std::string tagtype);
    ~AgvCalibration();

    /*
    boardSize: the number of inner corner of chessboard (row and coulum)
    boardSquareLength: side length of grid (m) 
    */
    void cameraUSB(int camera_index, 
                   cv::Size board_size,
                   string image_type);

    void cameraIP(string IP, 
                  cv::Size board_size,
                  string image_type);

    void imageToVector(int image_num, 
                       string image_path, 
                       string image_name, 
                       std::vector<cv::Mat> &images);

    void posToVector(string path_name, 
                     string file_name, 
                     int pos_num, 
                     std::vector<cv::Point3f> &objPosition);

    void intriCalib(std::vector<cv::Mat> chessBoardImgs,
                    cv::Size boardSize,
                    cv::Size2f boardSquareLenth,
                    cv::Mat &camera_matrix,
                    cv::Mat &camera_distortion);

    bool AprilTagDetection(TagShow &ts, 
                           cv::Mat &atimg, 
                           std::pair<float, float> &imageCenter, 
                           int index);

    bool extriCalib(std::vector<Mat> aprTagImages, 
                    std::vector<cv::Point3f> objPosition, 
                    cv::Mat intrixMatrix,
                    cv::Mat distCoff,
                    std::vector<Point2f> &tagCenters,
                    cv::Mat &rotation,
                    cv::Mat &transaction);
                    
    void cameraMatrixDoubleToFloat(cv::Mat &camera_matrix, 
                                   cv::Mat &camera_distortion);

    bool ImgPoints2World(std::vector<cv::Point2f> &pixels,
                         cv::Mat &cameraMatrix, 
                         cv::Mat &distCoeffs, 
                         cv::Mat &R_Cam2World, cv::Mat &t_Cam2World,
                         std::vector<cv::Point3f> &points,
                         double depth/* = -1.0*/, // Uints mm
                         int index/*=-1*/); // -1 for all.

    float reprojectError(std::vector<cv::Point3f> object_position, 
                         std::vector<cv::Point3f> reproject_position);

    void saveParameters(const string file_name,
                        const cv::Mat camera_matrix, 
                        const cv::Mat camera_distortion, 
                        const cv::Mat rotation, 
                        const cv::Mat translation);

    void loadParameters(const string file_name,
                        cv::Mat &camera_matrix, 
                        cv::Mat &camera_distortion, 
                        cv::Mat &rotation, 
                        cv::Mat &translation);

private:
    // cv::Mat intrixMatrix;
    // cv::Mat distCoff;
    // cv::Vec3d m_rvec, m_tvec;
    // std::vector<double> extrixParameters;
    double aprilTagSize;
    std::string tagType;
};



