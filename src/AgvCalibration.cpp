#include "AgvCalibration.h"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

/************************
 * discription:this function can convert images(jpg,bmp) to images vector(cv::Mat)
 * parameters: image_num:number of images
 *             image_path:folder path
 *             image_name:image name your saving
 *             images:output images vector
 * author:     chending
 * date:       2021/11/24
 *************************/
void AgvCalibration::imageToVector(int image_num,
                                   string image_path,
                                   string image_name,
                                   std::vector<cv::Mat> &images)
{
    for (int i = 1; i <= image_num; i++)
    {
        std::string file_path;
        file_path = image_path + image_name + to_string(i) + ".jpg";
        cv::Mat image = cv::imread(file_path, cv::IMREAD_COLOR);
        if (image.empty())
        {
            std::cout << "can not find picture" << std::endl;
            return;
        }
        else
        {
            images.push_back(image);
        }
    }
}

/************************
 * discription:this function can convert images(jpg,bmp) to images vector(cv::Mat)
 * parameters: image_num:number of images
 *             image_path:folder path
 *             image_name:image name your saving
 *             images:output images vector
 * author:     chending
 * date:       2021/11/24
 *************************/
void AgvCalibration::posToVector(string path_name,
                                 string file_name,
                                 int pos_num,
                                 std::vector<cv::Point3f> &objPosition)
{
    string file_input = path_name + file_name;
    std::ifstream fin;
    fin.open(file_input, std::ios::in);
    if (!fin.is_open())
    {
        std::cout << "pose file not found" << std::endl;
        return;
    }
    else
    {
        std::cout << "pose file found" << std::endl;
    }

    for (int i = 1; i <= pos_num; i++)
    {
        cv::Point3f one_point;
        double temp;

        // std::cout << "pose[" << i << "]:";
        fin >> temp;
        // std::cout << temp << "  " ;
        one_point.x = temp;
        fin >> temp;
        // std::cout << temp << "  " ;
        one_point.y = temp;
        fin >> temp;
        // std::cout << temp << "  "  << std::endl;
        one_point.z = temp;
        objPosition.push_back(one_point);
    }
}

void AgvCalibration::intriCalib(std::vector<cv::Mat> chess_board_imgs,
                                cv::Size boardSize,
                                cv::Size2f boardSquareLenth,
                                cv::Mat &camera_matrix,
                                cv::Mat &camera_distortion)
{
    std::cout << "the num of input chesseboard image: " << chess_board_imgs.size() << endl;
    std::vector<cv::Point2f> image_points;
    std::vector<std::vector<cv::Point2f>> image_points_buff;
    int imgcount = 0;
    cv::Size image_size;
    for (int i = 0; i < chess_board_imgs.size(); i++)
    {
        cv::Mat img = chess_board_imgs[i];
        if (i == 0)
        {
            image_size.width = img.cols;
            image_size.height = img.rows;
        }

        bool find = cv::findChessboardCorners(img, boardSize, image_points, cv::CALIB_CB_ADAPTIVE_THRESH | cv::CALIB_CB_NORMALIZE_IMAGE);
        if (find)
        {
            cv::cvtColor(img, img, cv::COLOR_BGR2GRAY);
            cv::cornerSubPix(img, image_points, cv::Size(5, 5), cv::Size(-1, -1), cv::TermCriteria(TermCriteria::EPS + TermCriteria::MAX_ITER, 50, DBL_EPSILON));
            // cv::find4QuadCornerSubpix(img, image_points, cv::Size(5, 5));
            image_points_buff.push_back(image_points);

            cv::Mat image_col;
            cv::cvtColor(img, image_col, COLOR_GRAY2BGR);
            cv::drawChessboardCorners(image_col, boardSize, image_points, true);
            cv::imwrite("../resource/result/chessboard_" + to_string(i) + ".jpg", image_col);
            imgcount++;
        }
        else
        {
            std::cout << "can not find chessboard corners: " << i << endl;
            continue;
        }
    }
    std::cout << "useful image size: " << imgcount << endl;

    int num_of_images = image_points_buff.size();
    std::vector<std::vector<cv::Point3f>> object_points;
    std::vector<int> point_counts;
    std::vector<cv::Mat> t_vec;
    std::vector<cv::Mat> r_vec;
    int i, j, k;
    for (k = 0; k < num_of_images; k++)
    {
        std::vector<cv::Point3f> temp_point_set;
        for (i = 0; i < boardSize.height; i++)
        {
            for (j = 0; j < boardSize.width; j++)
            {
                cv::Point3f real_point;
                real_point.x = j * boardSquareLenth.width;
                real_point.y = i * boardSquareLenth.height;
                real_point.z = 0;
                temp_point_set.push_back(real_point);
            }
        }
        object_points.push_back(temp_point_set);
    }
    // 初始化每幅图像上的角点数量
    for (int i = 0; i < num_of_images; i++)
    {
        point_counts.push_back(boardSize.width * boardSize.height);
    }

    cv::calibrateCamera(object_points, image_points_buff, image_size, camera_matrix, camera_distortion, r_vec, t_vec, 
                        0, TermCriteria(TermCriteria::MAX_ITER + TermCriteria::EPS, 100, DBL_EPSILON));

    // calculate re_project error
    double total_err = 0.0;
    double err = 0.0;

    std::vector<cv::Point2f> image_points2;
    for (int i = 0; i < num_of_images; i++)
    {
        std::vector<cv::Point3f> temp_point_set = object_points[i];

        // re_projected
        cv::projectPoints(temp_point_set, r_vec[i], t_vec[i], camera_matrix, camera_distortion, image_points2);

        std::vector<cv::Point2f> temp_image_points = image_points_buff[i];
        cv::Mat temp_image_points_Mat = cv::Mat(1, temp_image_points.size(), CV_32FC2);
        cv::Mat temp_image_points2_Mat = cv::Mat(1, image_points2.size(), CV_32FC2);
        for (int j = 0; j < temp_image_points.size(); j++)
        {
            temp_image_points_Mat.at<cv::Vec2f>(0, j) = cv::Vec2f(temp_image_points[j].x, temp_image_points[j].y);
            temp_image_points2_Mat.at<cv::Vec2f>(0, j) = cv::Vec2f(image_points2[j].x, image_points2[j].y);
        }
        err = cv::norm(temp_image_points_Mat, temp_image_points2_Mat, cv::NORM_L2);
        total_err += err /= point_counts[i];
        std::cout << i << "average reproject error:  " << err << "pixel" << std::endl;
    }

    std::cout << std::endl
              << "total rp error: " << total_err / num_of_images << "pixel" << std::endl;
}

void AgvCalibration::cameraMatrixDoubleToFloat(cv::Mat &camera_matrix,
                                               cv::Mat &camera_distortion)
{
    /*  double camera matrix covert to float camera matrix  */
    double camera_matrix_double00 = camera_matrix.at<double>(0, 0);
    float camera_matrix_float00 = (float)camera_matrix_double00;
    double camera_matrix_double01 = camera_matrix.at<double>(0, 1);
    float camera_matrix_float01 = (float)camera_matrix_double01;
    double camera_matrix_double02 = camera_matrix.at<double>(0, 2);
    float camera_matrix_float02 = (float)camera_matrix_double02;

    double camera_matrix_double10 = camera_matrix.at<double>(1, 0);
    float camera_matrix_float10 = (float)camera_matrix_double10;
    double camera_matrix_double11 = camera_matrix.at<double>(1, 1);
    float camera_matrix_float11 = (float)camera_matrix_double11;
    double camera_matrix_double12 = camera_matrix.at<double>(1, 2);
    float camera_matrix_float12 = (float)camera_matrix_double12;

    double camera_matrix_double20 = camera_matrix.at<double>(2, 0);
    float camera_matrix_float20 = (float)camera_matrix_double20;
    double camera_matrix_double21 = camera_matrix.at<double>(2, 1);
    float camera_matrix_float21 = (float)camera_matrix_double21;
    double camera_matrix_double22 = camera_matrix.at<double>(2, 2);
    float camera_matrix_float22 = (float)camera_matrix_double22;

    camera_matrix = (Mat_<float>(3, 3) << camera_matrix_float00, camera_matrix_float01, camera_matrix_float02,
                     camera_matrix_float10, camera_matrix_float11, camera_matrix_float12,
                     camera_matrix_float20, camera_matrix_float21, camera_matrix_float22);
    // std::cout << "camera_matrix: \n" << camera_matrix <<std::endl;
    // std::cout << "cameraMatrix: \n" << cameraMatrix <<std::endl;

    /*  double camera distortion covert to float camera distortion  */
    double camera_distortion_double00 = camera_distortion.at<double>(0, 0);
    float camera_distortion_float00 = (float)camera_distortion_double00;
    double camera_distortion_double01 = camera_distortion.at<double>(0, 1);
    float camera_distortion_float01 = (float)camera_distortion_double01;
    double camera_distortion_double02 = camera_distortion.at<double>(0, 2);
    float camera_distortion_float02 = (float)camera_distortion_double02;

    double camera_distortion_double03 = camera_distortion.at<double>(0, 3);
    float camera_distortion_float03 = (float)camera_distortion_double03;
    double camera_distortion_double04 = camera_distortion.at<double>(0, 4);
    float camera_distortion_float04 = (float)camera_distortion_double04;

    camera_distortion = (Mat_<float>(1, 5) << camera_distortion_float00, camera_distortion_float01, camera_distortion_float02,
                         camera_distortion_float03, camera_distortion_float04);

    // std::cout << "camera_distortion: \n" << camera_distortion <<std::endl;
    // std::cout << "distCoeff: \n" << distCoeff <<std::endl;
}

inline float GetSFromCalibration(cv::Point2f &pixel,
                                 cv::Mat &cameraMatrixInv,
                                 cv::Mat &R_Cam2World,
                                 cv::Mat &t_Cam2World)
{
    float pixels[3] = {pixel.x, pixel.y, 1.0};
    float r20 = R_Cam2World.at<float>(2, 0);
    float r21 = R_Cam2World.at<float>(2, 1);
    float r22 = R_Cam2World.at<float>(2, 2);
    float rVec[3] = {r20, r21, r22};
    float pz = t_Cam2World.at<float>(2);

    cv::Mat R2Vector(cv::Size(3, 1), CV_32FC1, rVec);
    cv::Mat pixelVector(cv::Size(1, 3), CV_32FC1, pixels);
    cv::Mat down = R2Vector * cameraMatrixInv * pixelVector;

    float upper = -1.0 * pz;

    return upper / down.at<float>(0);
}

bool AgvCalibration::ImgPoints2World(vector<cv::Point2f> &pixels,
                                     cv::Mat &cameraMatrix,
                                     cv::Mat &distCoeffs,
                                     cv::Mat &R_Cam2World, cv::Mat &t_Cam2World,
                                     vector<cv::Point3f> &points,
                                     double depth = -1.0, // Uints mm
                                     int index = -1)      // -1 for all.
{
    vector<cv::Point2f> undistorted;
    undistorted.resize(pixels.size());

    cv::undistortPoints(pixels, undistorted, cameraMatrix, distCoeffs /*, cv::noArray(), cameraMatrix*/);

    for (auto &pixel : undistorted)
    {
        pixel.x = pixel.x * cameraMatrix.at<float>(0, 0) + cameraMatrix.at<float>(0, 2);
        pixel.y = pixel.y * cameraMatrix.at<float>(1, 1) + cameraMatrix.at<float>(1, 2);

        cv::Point3f pixel3D{pixel.x, pixel.y, 1.0};
        cv::Mat cmInv, R_World2Cam;
        cv::invert(cameraMatrix, cmInv);
        cv::invert(R_Cam2World, R_World2Cam);
        cv::Mat point3DMat(pixel3D);

        if (depth < 0.0)
        {
            depth = GetSFromCalibration(pixel, cmInv, R_Cam2World, t_Cam2World);

            string info;
            if (isnan(depth))
            {
                info = "Calculated depth value is " + to_string(depth);
                // LogError(info.c_str());
                return false;
            }
        }

        cv::Mat point = (depth * R_Cam2World * cmInv * point3DMat) + t_Cam2World; // pixel to world
        cv::Point3f point3DWorld(point.at<float>(0), point.at<float>(1), point.at<float>(2));
        points.push_back(point3DWorld);
    }
    return false;
}

float AgvCalibration::reprojectError(std::vector<cv::Point3f> object_position,
                                     std::vector<cv::Point3f> reproject_position)
{
    float error_all = 0;
    if (object_position.size() != reproject_position.size())
    {
        std::cout << "size of object position is not match size of reproject position" << std::endl;
        return 0.0;
    }
    for (int i = 0; i < object_position.size(); i++)
    {
        float x_error = abs(object_position[i].x - reproject_position[i].x);
        float y_error = abs(object_position[i].y - reproject_position[i].y);
        float error = sqrt(x_error * x_error + y_error * y_error);
        // std::cout << "reproject error of point " << i << error <<std::endl;
        error_all += error;
    }

    float error_average = error_all / object_position.size();

    return (error_average);
}

void AgvCalibration::saveParameters(const string file_name,
                                    const cv::Mat camera_matrix,
                                    const cv::Mat camera_distortion,
                                    const cv::Mat R,
                                    const cv::Mat T)
{
    cv::FileStorage fs_write(file_name, cv::FileStorage::WRITE);

    fs_write << "camera_matrix" << camera_matrix;
    fs_write << "camera_distortion" << camera_distortion;
    fs_write << "R" << R;
    fs_write << "T" << T;

    fs_write.release();
}

void AgvCalibration::loadParameters(const string file_name,
                                    cv::Mat &camera_matrix,
                                    cv::Mat &camera_distortion,
                                    cv::Mat &R,
                                    cv::Mat &T)
{
    cv::FileStorage fs_read(file_name, cv::FileStorage::READ);

    fs_read["camera_matrix"] >> camera_matrix;
    fs_read["camera_distortion"] >> camera_distortion;
    fs_read["R"] >> R;
    fs_read["T"] >> T;

    fs_read.release();
}
