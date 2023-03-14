#include "calibration.h"

void Calibration::intrinsicCalibration(std::vector<cv::Mat> chess_board_imgs,
                                       cv::Size boardSize,
                                       cv::Size2f boardSquareLenth,
                                       cv::Mat &camera_matrix,
                                       cv::Mat &camera_distortion)
{
    std::cout << "the num of input chesseboard image: " << chess_board_imgs.size() << std::endl;
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
            cv::cornerSubPix(img, image_points, cv::Size(5, 5), cv::Size(-1, -1), cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::MAX_ITER, 50, DBL_EPSILON));
            // cv::find4QuadCornerSubpix(img, image_points, cv::Size(5, 5));
            image_points_buff.push_back(image_points);

            cv::Mat image_col;
            cv::cvtColor(img, image_col, cv::COLOR_GRAY2BGR);
            cv::drawChessboardCorners(image_col, boardSize, image_points, true);
            cv::imwrite("../resource/result/chessboard_" + std::to_string(i) + ".jpg", image_col);
            imgcount++;
        }
        else
        {
            std::cout << "can not find chessboard corners: " << i << std::endl;
            continue;
        }
    }
    std::cout << "useful image size: " << imgcount << std::endl;

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
                        0, cv::TermCriteria(cv::TermCriteria::MAX_ITER + cv::TermCriteria::EPS, 100, DBL_EPSILON));

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

// #include <iostream>

// #include <opencv2/opencv.hpp>

// uchar bilinearInter(double u, double v, const cv::Mat &image)
// {
//     if ((int)u >= 0 && (int)u <= image.cols - 2)
//     {
//         if ((int)v >= 0 && (int)v <= image.rows - 2)
//         {
//             uchar LT = image.at<uchar>((int)v, (int)u);
//             uchar LB = image.at<uchar>((int)v + 1, (int)u);
//             uchar RT = image.at<uchar>((int)v, (int)u + 1);
//             uchar RB = image.at<uchar>((int)v + 1, (int)u + 1);

//             double dLx, dTy, dRx, dBy;
//             dLx = u - (int)u;
//             dRx = (int)u - u + 1;
//             dTy = v - (int)v;
//             dBy = (int)v - v + 1;

//             uchar TP = dLx * RT + dRx * LT;
//             uchar TB = dLx * RB + dRx * LB;

//             return (TP * dBy + TB * dTy);
//         }
//     }
//     else
//         return 0;
// }

// cv::Mat &undistorted(const cv::Mat &init_image, const std::vector<double> intrinsic, const std::vector<double> distort)
// {
//     cv::Mat temp;
//     init_image.copyTo(temp);

//     if (intrinsic.size() + distort.size() < 8)
//     {
//         std::cout << "error" << std::endl;
//     }

//     // double intrinisic[4] = { 458.654,457.296,367.215,248.375 };                                       //double intrinsic[4]={fx,fy,cx,cy};
//     // double disort[4] = { -0.28340811,0.07395907,0.00019359,1.76187114e-05 };                          //double disort[4]={k1,k2,p1,p2};

//     cv::Mat distorted_image = cv::Mat(init_image.rows, init_image.cols, init_image.type());
//     cv::Mat nearest_inter_image = cv::Mat(init_image.rows, init_image.cols, init_image.type());
//     for (int v = 0; v < init_image.rows; v++)
//     {
//         for (int u = 0; u < init_image.cols; u++)
//         {
//             /* 像素平面 >> 归一化平面 */
//             double X_1 = (u - intrinsic.at(2)) / intrinsic.at(0);
//             double Y_1 = (v - intrinsic.at(3)) / intrinsic.at(1);

//             /* 归一化平面去畸变 */
//             double R = sqrt(X_1 * X_1 + Y_1 * Y_1);
//             double X_undistorted = X_1 * (1 + distort.at(0) * R * R + distort.at(1) * R * R * R * R) + 2 * distort.at(2) * X_1 * Y_1 + distort.at(3) * (R * R + 2 * X_1 * X_1);
//             double Y_undistorted = Y_1 * (1 + distort.at(0) * R * R + distort.at(1) * R * R * R * R) + 2 * distort.at(2) * (R * R + 2 * Y_1 * Y_1) + 2 * distort.at(3) * X_1 * Y_1;

//             /* 归一化平面 >> 像素平面(去畸变后) */
//             double U_undistorted = intrinsic.at(0) * X_undistorted + intrinsic.at(2);
//             double V_undistorted = intrinsic.at(1) * Y_undistorted + intrinsic.at(3);

//             /* 像素平面优化（插值填充） */
//             distorted_image.at<uchar>(v, u) = bilinearInter(U_undistorted, V_undistorted, init_image); // 双线性插值计算去畸变的像素值；

//             /* 临近插值 */
//             if (U_undistorted >= 0 && V_undistorted >= 0 && U_undistorted < init_image.cols && V_undistorted < init_image.rows)
//                 nearest_inter_image.at<uchar>(v, u) = init_image.at<uchar>((int)V_undistorted, (int)U_undistorted);
//             else
//                 nearest_inter_image.at<uchar>(v, u) = 0;
//         }
//     }

//     cv::imshow("nearestInterpolation", nearest_inter_image);
//     cv::imshow("bilinearInterpolation", distorted_image);
//     cv::imshow("initialImage", init_image);
//     cv::imwrite("nearestInterpolation.png", nearest_inter_image);
//     cv::imwrite("bilinearInterpolation.png", distorted_image);
//     cv::waitKey(0);
//     return distorted_image;
// }

// int main()
// {
//     cv::Mat image = cv::imread("../chess_8.jpg");
//     std::vector<double>intrinsic({837.2869180323291, 841.8236237071096, 556.3932298660817, 387.4635354570347});
//     std::vector<double>distort({-0.4480429340370656, 0.2874356149466112, 0.002370861170220438, 0.003536092926783489});
//     undistorted(image, intrinsic, distort);

//     return 0;
// }

// #include <opencv2/opencv.hpp>
// #include <string>
// #include <math.h>

// using namespace std;

// uchar Bilinear_Inter(double u, double v, const cv::Mat &image)
// {
//     if ((int)u >= 0 && (int)u <= image.cols - 2)
//         if ((int)v >= 0 && (int)v <= image.rows - 2)
//         {
//             uchar
//                 LT = image.at<uchar>((int)v, (int)u),
//                 LB = image.at<uchar>((int)v + 1, (int)u),
//                 RT = image.at<uchar>((int)v, (int)u + 1),
//                 RB = image.at<uchar>((int)v + 1, (int)u + 1);

//             double dLx, dTy, dRx, dBy;
//             dLx = u - (int)u;
//             dRx = (int)u - u + 1;
//             dTy = v - (int)v;
//             dBy = (int)v - v + 1;
//             uchar TP = dLx * RT + dRx * LT;

//             uchar TB = dLx * RB + dRx * LB;
//             return TP * dBy + TB * dTy;
//         }
//         else
//             return 0;
// }

// string image_file = "../chess_8.jpg"; // 请确保路径正确

// int main(int argc, char **argv)
// {

//     // 本程序需要你自己实现去畸变部分的代码。尽管我们可以调用OpenCV的去畸变，但自己实现一遍有助于理解。
//     // 畸变参数
//     double k1 = -0.4480429340370656, k2 = 0.2874356149466112, p1 = 0.002370861170220438, p2 = 0.003536092926783489, k3 = -0.1202712527495983; //-0.4480429340370656, 0.2874356149466112, 0.002370861170220438, 0.003536092926783489, -0.1202712527495983
//     // 内参
//     double fx = 837.2869180323291, fy = 841.8236237071096, cx = 556.3932298660817, cy = 841.8236237071096;

//     cv::Mat image = cv::imread(image_file, 0); // 图像是灰度图，CV_8UC1
//     int rows = image.rows, cols = image.cols;
//     cv::Mat image_undistort = cv::Mat(rows, cols, CV_8UC1); // 去畸变以后的图

//     // 计算去畸变后图像的内容
//     for (int v = 0; v < rows; v++)
//     {
//         for (int u = 0; u < cols; u++)
//         {

//             double u_distorted = 0, v_distorted = 0;
//             // TODO 按照公式，计算点(u,v)对应到畸变图像中的坐标(u_distorted, v_distorted) (~6 lines)
//             // start your code here
//             // image_undistort中含有非畸变的图像坐标
//             // 将image_undistort的坐标通过内参转换到归一化坐标系下，此时得到的归一化坐标是对的
//             // 将得到的归一化坐标系进行畸变处理
//             // 将畸变处理后的坐标通过内参转换为图像坐标系下的坐标
//             // 这样就相当于是在非畸变图像的图像坐标和畸变图像的图像坐标之间建立了一个对应关系
//             // 相当于是非畸变图像坐标在畸变图像中找到了映射
//             // 对畸变图像进行遍历之后，然后赋值（一般需要线性插值，因为畸变后图像的坐标不一定是整数的），即可得到矫正之后的图像
//             double x1, y1, x2, y2;
//             x1 = (u - cx) / fx;
//             y1 = (v - cy) / fy;
//             double r2;
//             r2 = pow(x1, 2) + pow(y1, 2);
//             x2 = x1 * (1 + k1 * r2 + k2 * pow(r2, 2) + k3 * pow(r2, 3)) + 2 * p1 * x1 * y1 + p2 * (r2 + 2 * x1 * x1);
//             y2 = y1 * (1 + k1 * r2 + k2 * pow(r2, 2) + k3 * pow(r2, 3)) + p1 * (r2 + 2 * y1 * y1) + 2 * p2 * x1 * y1;

//             u_distorted = fx * x2 + cx;
//             v_distorted = fy * y2 + cy;

//             // end your code here
//             image_undistort.at<uchar>(v, u) = Bilinear_Inter(u_distorted, v_distorted, image);

            // 赋值 (最近邻插值)
            // if (u_distorted >= 0 && v_distorted >= 0 && u_distorted < cols && v_distorted < rows)
            // {
            //     image_undistort.at<uchar>(v, u) = image.at<uchar>((int)v_distorted, (int)u_distorted);
            // }
            // else
            // {
            //     image_undistort.at<uchar>(v, u) = 0;
            // }
//         }
//     }
//     // 画图去畸变后图像
//     cv::imshow("image undistorted", image_undistort);
//     cv::waitKey();

//     return 0;
// }

// #include <opencv2/opencv.hpp>
// #include <string>

// using namespace std;

// string image_file = "../chess_8.jpg"; // 请确保路径正确

// int main(int argc, char **argv)
// {

//     // 本程序实现去畸变部分的代码。尽管我们可以调用OpenCV的去畸变，但自己实现一遍有助于理解。
//     // 畸变参数
//     double k1 = -0.4480429340370656, k2 = 0.2874356149466112, p1 = 0.002370861170220438, p2 = 0.003536092926783489, k3 = -0.1202712527495983;
//     // 内参
//     double fx = 837.2869180323291, fy = 841.8236237071096, cx = 556.3932298660817, cy = 841.8236237071096;

//     cv::Mat image = cv::imread(image_file, 0); // 图像是灰度图，CV_8UC1
//     int rows = image.rows, cols = image.cols;
//     cv::Mat image_undistort = cv::Mat(rows, cols, CV_8UC1); // 去畸变以后的图

//     // 计算去畸变后图像的内容
//     for (int v = 0; v < rows; v++)
//     {
//         for (int u = 0; u < cols; u++)
//         {
//             // 按照公式，计算点(u,v)对应到畸变图像中的坐标(u_distorted, v_distorted)
//             double x = (u - cx) / fx, y = (v - cy) / fy;
//             double r = sqrt(x * x + y * y);

//             // double x_undistorted, y_undistorted;
//             // x_distorted = x_undistorted * (1 + k1 * r * r + k2 * r * r * r * r) + 2 * p1 * x_undistorted * y_undistorted + p2 * (r * r + 2 * x_undistorted * x_undistorted);
//             // y_distorted = y_undistorted * (1 + k1 * r * r + k2 * r * r * r * r) + p1 * (r * r + 2 * y_undistorted * y_undistorted) + 2 * p2 * x_undistorted * y_undistorted;
//             // double u_undistorted, v_undistorted;
//             // u_undistorted = fx * x_undistorted + cx;
//             // v_undistorted = fy * y_undistorted + cy;

//             double x_distorted = x * (1 + k1 * r * r + k2 * r * r * r * r) + 2 * p1 * x * y + p2 * (r * r + 2 * x * x);
//             double y_distorted = y * (1 + k1 * r * r + k2 * r * r * r * r) + p1 * (r * r + 2 * y * y) + 2 * p2 * x * y;
//             double u_distorted = fx * x_distorted + cx;
//             double v_distorted = fy * y_distorted + cy;

//             // 赋值 (最近邻插值)
//             if (u_distorted >= 0 && v_distorted >= 0 && u_distorted < cols && v_distorted < rows)
//             {
//                 image_undistort.at<uchar>(v, u) = image.at<uchar>((int)v_distorted, (int)u_distorted);
//             }
//             else
//             {
//                 image_undistort.at<uchar>(v, u) = 0;
//             }
//         }
//     }

//     // 画图去畸变后图像
//     cv::imshow("distorted", image);
//     cv::imshow("undistorted", image_undistort);
//     cv::waitKey();
//     return 0;
// }
