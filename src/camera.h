/*
 * @Author: chending2020529 chending529@gmail.com
 * @Date: 2023-03-13 20:56:00
 * @LastEditors: chending2020529 chending529@gmail.com
 * @LastEditTime: 2023-03-13 21:57:39
 * @FilePath: /LearningOpenCV/src/camera.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#pragma once

#include <iostream>

#include <opencv2/opencv.hpp>

class Cameras
{
public:
    Cameras() = default;
    ~Cameras() = default;
    Cameras(const Cameras &cam);
    Cameras &operator=(Cameras &cam);

public:
    void cameraUSB(int camera_index, cv::Mat &image);
    void cameraRTSP(std::string IP, cv::Mat &image);
};