// /*
//  * @Author: chending chending2@seer-group.com
//  * @Date: 2023-03-13 11:05:34
//  * @LastEditors: chending chending2@seer-group.com
//  * @LastEditTime: 2023-03-13 12:32:07
//  * @FilePath: /LearningOpenCV/src/calibration.h
//  * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
//  */
// #pragma once

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

//             // 赋值 (最近邻插值)
//             // if (u_distorted >= 0 && v_distorted >= 0 && u_distorted < cols && v_distorted < rows)
//             // {
//             //     image_undistort.at<uchar>(v, u) = image.at<uchar>((int)v_distorted, (int)u_distorted);
//             // }
//             // else
//             // {
//             //     image_undistort.at<uchar>(v, u) = 0;
//             // }
//         }
//     }
//     // 画图去畸变后图像
//     cv::imshow("image undistorted", image_undistort);
//     cv::waitKey();

//     return 0;
// }

#include <opencv2/opencv.hpp>
#include <string>

using namespace std;

string image_file = "../chess_8.jpg"; // 请确保路径正确

int main(int argc, char **argv)
{

    // 本程序实现去畸变部分的代码。尽管我们可以调用OpenCV的去畸变，但自己实现一遍有助于理解。
    // 畸变参数
    double k1 = -0.4480429340370656, k2 = 0.2874356149466112, p1 = 0.002370861170220438, p2 = 0.003536092926783489, k3 = -0.1202712527495983;
    // 内参
    double fx = 837.2869180323291, fy = 841.8236237071096, cx = 556.3932298660817, cy = 841.8236237071096;

    cv::Mat image = cv::imread(image_file, 0); // 图像是灰度图，CV_8UC1
    int rows = image.rows, cols = image.cols;
    cv::Mat image_undistort = cv::Mat(rows, cols, CV_8UC1); // 去畸变以后的图

    // 计算去畸变后图像的内容
    for (int v = 0; v < rows; v++)
    {
        for (int u = 0; u < cols; u++)
        {
            // 按照公式，计算点(u,v)对应到畸变图像中的坐标(u_distorted, v_distorted)
            double x = (u - cx) / fx, y = (v - cy) / fy;
            double r = sqrt(x * x + y * y);

            // double x_undistorted, y_undistorted;
            // x_distorted = x_undistorted * (1 + k1 * r * r + k2 * r * r * r * r) + 2 * p1 * x_undistorted * y_undistorted + p2 * (r * r + 2 * x_undistorted * x_undistorted);
            // y_distorted = y_undistorted * (1 + k1 * r * r + k2 * r * r * r * r) + p1 * (r * r + 2 * y_undistorted * y_undistorted) + 2 * p2 * x_undistorted * y_undistorted;
            // double u_undistorted, v_undistorted;
            // u_undistorted = fx * x_undistorted + cx;
            // v_undistorted = fy * y_undistorted + cy;

            double x_distorted = x * (1 + k1 * r * r + k2 * r * r * r * r) + 2 * p1 * x * y + p2 * (r * r + 2 * x * x);
            double y_distorted = y * (1 + k1 * r * r + k2 * r * r * r * r) + p1 * (r * r + 2 * y * y) + 2 * p2 * x * y;
            double u_distorted = fx * x_distorted + cx;
            double v_distorted = fy * y_distorted + cy;

            // 赋值 (最近邻插值)
            if (u_distorted >= 0 && v_distorted >= 0 && u_distorted < cols && v_distorted < rows)
            {
                image_undistort.at<uchar>(v, u) = image.at<uchar>((int)v_distorted, (int)u_distorted);
            }
            else
            {
                image_undistort.at<uchar>(v, u) = 0;
            }
        }
    }

    // 画图去畸变后图像
    cv::imshow("distorted", image);
    cv::imshow("undistorted", image_undistort);
    cv::waitKey();
    return 0;
}
