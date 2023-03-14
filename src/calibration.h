/*
 * @Author: chending2020529 chending529@gmail.com
 * @Date: 2023-03-13 21:55:47
 * @LastEditors: chending2020529 chending529@gmail.com
 * @LastEditTime: 2023-03-14 08:28:51
 * @FilePath: /LearningOpenCV/src/calibration.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#pragma once
#include <vector>
#include <opencv2/opencv.hpp>

class Calibration
{
public:
    Calibration() = default;
    ~Calibration() = default;
    Calibration(const Calibration &calib);
    Calibration &operator=(const Calibration &calib);

public:
    void intrinsicCalibration(std::vector<cv::Mat> chessBoardImgs,
                              cv::Size boardSize,
                              cv::Size2f boardSquareLenth,
                              cv::Mat &camera_matrix,
                              cv::Mat &camera_distortion);
};