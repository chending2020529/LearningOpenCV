/*
 * @Author: chending2020529 chending529@gmail.com
 * @Date: 2023-03-13 20:41:27
 * @LastEditors: chending2020529 chending529@gmail.com
 * @LastEditTime: 2023-03-13 23:11:04
 * @FilePath: /LearningOpenCV/src/AgvCalibration.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>

using namespace std;
using namespace cv;

class AgvCalibration
{
public:
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



